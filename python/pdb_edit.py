'''
Useful manipulations on PDB files
'''

# Python imports
import copy
import os
import re

# our imports
import ample_util
import pdb_model


three2one = {
    'ALA' : 'A',    
    'ARG' : 'R',    
    'ASN' : 'N',    
    'ASP' : 'D',    
    'CYS' : 'C',    
    'GLU' : 'E',    
    'GLN' : 'Q',    
    'GLY' : 'G',    
    'HIS' : 'H',    
    'ILE' : 'I',    
    'LEU' : 'L',    
    'LYS' : 'K',    
    'MET' : 'M',    
    'PHE' : 'F',    
    'PRO' : 'P',    
    'SER' : 'S',    
    'THR' : 'T',    
    'TRP' : 'W',    
    'TYR' : 'Y',   
    'VAL' : 'V'
}

# http://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python
#aaDict.update( dict((v, k) for (k, v) in aaDict.items()) )
one2three =  dict((v, k) for (k, v) in three2one.items()) 

class PDBEdit(object):
    """Class for editing PDBs
    
    """
    
#     def backbone(self, inpath=None, outpath=None ):
#         """Only output backbone atoms.
#         """        
#         
#         atom_names = [ 'N', 'CA', 'C', 'O', 'CB' ]
# 
#         #   print 'Found ',each_file
#         pdb_in = open( inpath, "r" )
#         pdb_out = open( outpath, "w" )    
# 
#         for pdbline in pdb_in:
#             pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
#             pdb_result = pdb_pattern.match(pdbline)
#     
#             if pdb_result:
#                 pdb_result2 = re.split(pdb_pattern, pdbline)
#                 if pdb_result2[3] != '':
#                     if pdb_result2[2] not in atom_names:
#                         continue
#             
#             # Write out everything else
#             pdb_out.write(pdbline)
#         
#         #End for
#         pdb_out.close()
#         pdb_in.close()
#         
#         return
    
    def backbone(self, inpath=None, outpath=None ):
        """Only output backbone atoms.
        """        
        
        logfile = outpath+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( inpath, outpath ).split()
        
        # Build up stdin
        stdin='lvatom "N,CA,C,O,CB[N,C,O]"'
        #stdin='lvatom "N,CA,C,O[N,C,O]"'
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
        
        if retcode == 0:
            # remove temporary files
            os.unlink(logfile)
        else:
            raise RuntimeError,"Error stripping PDB to backbone atoms"
        
        return
    
    def calpha_only( self, inpdb, outpdb ):
        """Strip PDB to c-alphas only"""
        
        
        logfile = outpdb+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( inpdb, outpdb ).split()
        
        # Build up stdin
        stdin='lvatom "CA[C]:*"'
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
        
        if retcode == 0:
            # remove temporary files
            os.unlink(logfile)
        else:
            raise RuntimeError,"Error stripping PDB to c-alpha atoms"
            
        return
    
    
    def extract_chain( self, inpdb, outpdb, chainID=None, newChainID=None, cAlphaOnly=False ):
        """Extract chainID from inpdb and renumner.
        If cAlphaOnly is set, strip down to c-alpha atoms
        """
        
        
        logfile = outpdb+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( inpdb, outpdb ).split()
        
        # Build up stdin
        stdin="lvchain {0}\n".format( chainID )
        if newChainID:
            stdin += "renchain {0} {1}\n".format( chainID, newChainID )
        if cAlphaOnly:
            stdin += 'lvatom "CA[C]:*"\n'
        stdin += "sernum\n"
        
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
        
        if retcode == 0:
            # remove temporary files
            os.unlink(logfile)
        else:
            raise RuntimeError,"Error extracting chain {0}".format( chainID )
            
        return

    def extract_model( self, inpdb, outpdb, modelID=None ):
        """Extract modelID from inpdb into outpdb"""
        
        assert modelID
        
        logfile = outpdb+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( inpdb, outpdb ).split()
        
        # Build up stdin
        stdin="lvmodel /{0}\n".format( modelID )
        #stdin += "sernum\n"
        
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
        
        if retcode != 0:
            raise RuntimeError,"Problem extracting model with cmd: {0}".format
    
        # remove temporary files
        os.unlink(logfile)
            
        return
    

    def get_resseq_map( self, nativePdb, modelPdb ):
        """Return a ResSeqMap mapping the index of a residue in the model to the corresponding residue in the native.
        Only works if 1 chain in either file and with standard residues
        """
        
        def _get_indices( pdb ):
            """Get sequence as string of 1AA
            get list of matching resSeq
            """
            
            
            print "GETTING INDICES ",pdb
            
            sequence = ""
            resSeq = []
            
            atomTypes = [] # For checking we have all required atom types
            backbone = [ 'N', 'CA', 'C', 'O','CB' ]
            
            backboneMask = []
            cAlphaMask = []
            
            chain=None
            readingResSeq=None
            readingResName=None
            for line in open( pdb ):
                
                if line.startswith("MODEL"):
                    raise RuntimeError,"FOUND MULTI_MODEL FILE!"
                
                if line.startswith("ATOM"):
                    
                    atom = pdb_model.PdbAtom( line )
                    
                    if not chain:
                        chain = atom.chainID
                    
                    if atom.chainID != chain:
                        raise RuntimeError," FOUND ADDITIONAL CHAIN"
                        break
                    
                    if atom.name not in atomTypes:
                        atomTypes.append( atom.name.strip() )
                        
                    # First atom in first residue
                    if readingResSeq == None:
                        readingResSeq = atom.resSeq
                        readingResName = atom.resName
                        continue
                    
                    if readingResSeq != atom.resSeq:
                        # Adding a new residue
                        
                        

                        
                        # Add the atom we've just finished reading
                        sequence += three2one[ readingResName ]
                        resSeq.append( readingResSeq )
                        
                        got=False
                        if 'CA' not in atomTypes:
                            cAlphaMask.append( True )
                            got=True
                        else:
                            cAlphaMask.append( False )
                            
                        
                        if not got: # If we haven't got CA we don't need to check
                            for at in backbone: # If we need to mask this residue for backbone atoms
                                if at not in atomTypes:
                                    got=True
                                    break
                                    #s = "Atom type {0} is not present in atom types for residue {1} - only got atomTypes: {2}".format( at, readingResSeq, atomTypes )
                                    #print s
                                    
                        if got:
                            backboneMask.append( True )
                        else:
                            backboneMask.append( False )
                            
                        
                        readingResSeq = atom.resSeq
                        readingResName = atom.resName
                        atomTypes = []
                        
            return ( sequence, resSeq, cAlphaMask, backboneMask )
      
        native_seq, native_idx = _get_indices( nativePdb )
        model_seq, model_idx = _get_indices( modelPdb )
        
        # The window of AA we used to check for a match    
        PROBE_LEN = 10
        
        if len(native_seq) < 20 or len(model_seq) < 20:
            raise RuntimeError,"Very short sequences - this will not work!"
        
        # MAXINSET is the max number of AA into the sequence that we will go searching for a match - i.e. if more
        # then MAXINSET AA are non-matching, we won't find the match 
        #MAXINSET=30 if len( model_seq ) > 30 else len( model_seq ) - ( PROBE_LEN + 2)
        if len( model_seq ) > 30:
            MAXINSET=30
        else:
            MAXINSET = len( model_seq ) - ( PROBE_LEN + 2)

        got=False
        for model_i in range( MAXINSET ):
            probe = model_seq[ model_i : model_i+PROBE_LEN-1 ]
            for native_i in range( MAXINSET ):
                if native_seq[ native_i:native_i+PROBE_LEN-1 ] == probe:
                    got=True
                    break
            
            if got:
                #print "GOT MODEL MATCH AT i,j ",model_i,native_i
                break
        
        # Now we know where they start we can sort out the indicies
        # map goes from the model -> native. For any in model that are not in native we set them to None
        resMap = residueSequenceMap()
        resMap.modelResSeq = model_idx
        
        for i in range( len( model_seq ) ):
            
            if i < model_i:
                # These are residues that are present in the model but not in the native
                resMap.nativeResSeq.append( None )
                continue
            
            pos = i - model_i + native_i
            if pos >= len( native_idx ):
                resMap.nativeResSeq.append(  None  )
            else:
                resMap.nativeResSeq.append(  native_idx[ pos ]  )
                
        
        if resMap.nativeResSeq != resMap.modelResSeq:
            raise RuntimeError, "Mismatching maps: {0}".format( resMap ) 
            
        return resMap

    def keep_matching( self, refpdb=None, targetpdb=None, outpdb=None, resSeqMap=None ):
        """Only keep those atoms in targetpdb that are in refpdb and write the result to outpdb.
        We also take care of renaming any chains.
        """
        
        assert refpdb and targetpdb and outpdb and resSeqMap
    
        # Paranoid check
        refinfo = self.get_info( refpdb )
        targetinfo = self.get_info( targetpdb )
        if len(refinfo.models) > 1 or len(targetinfo.models) > 1:
            raise RuntimeError, "PDBS contain more than 1 model!"
        
        if refinfo.models[0].chains != targetinfo.models[0].chains:
            raise RuntimeError, "Different numbers/names of chains {0}->{1} between {2} and {3}!".format( refinfo.models[0].chains,
                                                                                                        targetinfo.models[0].chains,
                                                                                                        refpdb,
                                                                                                        targetpdb
                                                                                                        )
        # Now we do our keep matching    
        tmp1 = ample_util.tmpFileName()+".pdb" # pdbcur insists names have a .pdb suffix
        
        self._keep_matching( refpdb, targetpdb, tmp1, resSeqMap=resSeqMap )
        
        # now renumber with pdbcur
        logfile = tmp1+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( tmp1, outpdb ).split()
        stdint="""sernum
    """
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdint)
        
        if retcode == 0:
            # remove temporary files
            os.unlink(tmp1)
            os.unlink(logfile)
        
        return retcode

    def _keep_matching( self, refpdb=None, targetpdb=None, outpdb=None, resSeqMap=None ):
        """Create a new pdb file that only contains that atoms in targetpdb that are
        also in refpdb. It only considers ATOM lines and discards HETATM lines in the target.
        
        Args:
        refpdb: path to pdb that contains the minimal set of atoms we want to keep
        targetpdb: path to the pdb that will be stripped of non-matching atoms
        outpdb: output path for the stripped pdb
        """
        
        assert refpdb and targetpdb and outpdb and resSeqMap
        
        def _output_residue( refResidues, targetAtomList, resSeqMap, outfh ):
            """Output a single residue only outputting matching atoms, shuffling the atom order and changing the resSeq num""" 
            
            # Get the matching list of atoms
            targetResSeq = targetAtomList[0].resSeq
            
            refResSeq = resSeqMap.native2model( targetResSeq )
            
            # Get the atomlist for the reference
            for ( rid, alist ) in refResidues:
                if rid == refResSeq:
                    refAtomList = alist
                    break
                
            # Get ordered list of the ref atom names for this residue
            rnames = [ x.name for x in refAtomList ]
            #print "got rnames ",rnames
            #print "got anames ", [ x.name for x in targetAtomList ]
            
            if len( refAtomList ) > len( targetAtomList ):
                s = "Cannot keep matching as refAtomList is > targetAtomList for residue {0}\nRef: {1}\nTrg: {2}".format( targetResSeq,
                                                                                                                          rnames,
                                                                                                                           [ x.name for x in targetAtomList ]
                                                                                                                           )
                raise RuntimeError, s

            # Remove any not matching in the target
            alist = []
            for atom in targetAtomList:
                if atom.name in rnames:
                    alist.append( atom )
            
            # List now only contains matching atoms
            targetAtomList = alist
            #print "tnames ",[ x.name for x in targetAtomList ]

            # Now just have matching so output in the correct order
            for refname in rnames:
                for i,atom in enumerate( targetAtomList ):
                    if atom.name == refname:
                        # Found the matching atom
                        
                        # Change resSeq and write out
                        atom.resSeq = refResSeq
                        outfh.write( atom.toLine()+"\n" )
                        
                        # now delete both this atom and the line
                        targetAtomList.pop(i)
                        
                        # jump out of inner loop
                        break   
            return
        
        # Go through refpdb and find which refResidues are present
        # ordered list of tuples - ( resSeq, [ list_of_atoms_for_that_residue ] )
        refResidues = []
        targetResSeq = []
        
        last = None
        chain = -1
        for line in open(refpdb, 'r'):
            
            if line.startswith("MODEL"):
                raise RuntimeError, "Multi-model file!"
            
            if line.startswith("TER"):
                break
            
            if line.startswith("ATOM"):
                a = pdb_model.PdbAtom( line )
                
                # First atom/chain
                if chain == -1:
                    chain = a.chainID
                
                if a.chainID != chain:
                    raise RuntimeError, "ENCOUNTERED ANOTHER CHAIN! {0}".format( line )

                if a.resSeq != last:
                    last = a.resSeq
                    
                    # Add the corresponding resSeq in the target
                    targetResSeq.append( resSeqMap.model2native( a.resSeq ) )
                    refResidues.append( ( a.resSeq, [ a ] ) )
                else:
                    refResidues[ -1 ][ 1 ].append( a )
                    
        # Now read in target pdb and output everything bar the atoms in this file that
        # don't match those in the refpdb
        t = open(targetpdb,'r')
        out = open(outpdb,'w')
        
        chain=None # The chain we're reading
        residue=None # the residue we're reading
        targetResidue = []
        
        for line in t:
            
            if line.startswith("MODEL"):
                raise RuntimeError, "Multi-model file!"

            if line.startswith("ANISOU"):
                raise RuntimeError, "I cannot cope with ANISOU! {0}".format(line)
            
            # Stop at TER
            if line.startswith("TER"):
                _output_residue( refResidues, targetResidue, resSeqMap, out )
                # we write out our own TER
                out.write("TER\n")
                continue
            
            if line.startswith("ATOM"):
                
                atom = pdb_model.PdbAtom( line )

                # First atom/chain
                if chain == None:
                    chain = a.chainID
                
                if atom.chainID != chain:
                    raise RuntimeError, "ENCOUNTERED ANOTHER CHAIN! {0}".format( line )
                
                if atom.resSeq in targetResSeq:
                    
                    # If this is the first one add the empty tuple and reset residue
                    if atom.resSeq != residue:
                        if residue != None: # Dont' write out owt for first atom
                            _output_residue( refResidues, targetResidue, resSeqMap, out )
                        targetResidue = []
                        residue = atom.resSeq
                    
                    # If not first keep adding
                    targetResidue.append( atom )
                    
                    # We don't write these out as we write them with _output_residue
                    continue
                    
                else:
                    # discard this line as not a matching atom
                    continue
            
            # For time being exclude all HETATM lines
            elif line.startswith("HETATM"):
                continue
            #Endif line.startswith("ATOM")
            
            # Output everything else
            out.write(line)
            
        # End reading loop
        
        t.close()
        out.close()
        
        return

# OLD STUFF FOR MULTIPLE CHAINS AND WITHOUT THE residueSequenceMap
#     def _keep_matching( self, refpdb=None, targetpdb=None, outpdb=None ):
#         """Create a new pdb file that only contains that atoms in targetpdb that are
#         also in refpdb. It only considers ATOM lines and discards HETATM lines in the target.
#         
#         Args:
#         refpdb: path to pdb that contains the minimal set of atoms we want to keep
#         targetpdb: path to the pdb that will be stripped of non-matching atoms
#         outpdb: output path for the stripped pdb
#         """
#     
#         assert refpdb and targetpdb and outpdb
#         
#         def _write_matching_residues( chain, refResidues, target_residues, outfh ):
#             
#             #print "got target_residues: {0}".format(target_residues)
#             
#             # Loop over each residue in turn
#             for idx, atoms_and_lines  in sorted( target_residues[ chain ].items() ):
#                 
#                 # Get ordered list of the ref atom names for this residue
#                 rnames = [ x.name for x in refResidues[ chain ][ idx ] ]
#                 
#                 #print "rnames ",rnames
#                 
#                 # Remove any not matching
#                 atoms = []
#                 atom_lines = []
#                 for i, a in enumerate( atoms_and_lines[0] ):
#                     if a.name in rnames:
#                         atoms.append( atoms_and_lines[0][i] )
#                         atom_lines.append( atoms_and_lines[1][i] )
#                 
#                 
#                 # Now just have matching so output in the correct order
#                 for refname in rnames:
#                     for i, atom in enumerate( atoms ):
#                         if atom.name == refname:
#                             # Found the matching atom so write out the corresponding line
#                             outfh.write( atom_lines[i] )
#                             # now delete both this atom and the line
#                             atoms.pop(i)
#                             atom_lines.pop(i)
#                             # jump out of inner loop
#                             break
#                         
#             # We delete the chain we've written out so that we don't write it out again at the
#             # end by mistake
#             del refResidues[ chainIdx ]
#             del target_residues[ chainIdx ]
#             return
#     
#         # Go through refpdb and find which refResidues are present
#         f = open(refpdb, 'r')
#         
#         # map of resSeq to list of PdbAtom objects for the reference residues
#         refResidues = {}
#         
#         last = None
#         chain = -1
#         chainIdx=-1 # For the time being we key by the chain index so we can deal with 
#                     # proteins that have different chain IDs
#         for line in f:
#             if line.startswith("MODEL"):
#                 raise RuntimeError, "Multi-model file!"
#             
#             if line.startswith("ATOM"):
#                 a = pdb_model.PdbAtom( line )
#                 
#                 if a.chainID != chain:
#                     chain = a.chainID
#                     chainIdx+=1
#                     if chainIdx in refResidues:
#                         raise RuntimeError, "ENCOUNTERED CHAIN AGAIN! {0}".format( line )
#                     refResidues[ chainIdx ] = {}
#                 
#                 if a.resSeq != last:
#                     #if a.resSeq in refResidues:
#                     #    raise RuntimeError,"Multiple chains in pdb - found residue #: {0} again.".format(a.resSeq)
#                     last = a.resSeq
#                     #refResidues[ last ] = [ a ]
#                     refResidues[ chainIdx ][ last ] = [ a ]
#                 else:
#                     #refResidues[ last ].append( a )
#                     refResidues[ chainIdx ][ last ].append( a )
#                     
#         f.close()
#         
#         #print "got refResidues: {0}".format(refResidues)
#         
#         # Now read in target pdb and output everything bar the atoms in this file that
#         # don't match those in the refpdb
#         t = open(targetpdb,'r')
#         out = open(outpdb,'w')
#         
#         reading=-1 # The residue we are reading - set to -1 when we are not reading
#         chain=-1 # The chain we're reading
#         chainIdx=-1 # see above
#         
#         target_residues = {} # dict mapping residue index to a a tuple of (atoms, lines), where atoms is a list of the atom
#         # objects and lines is a list of the lines used to create the atom objects
#         
#         for line in t:
#             
#             if line.startswith("MODEL"):
#                 raise RuntimeError, "Multi-model file!"
# 
#             if line.startswith("ANISOU"):
#                 raise RuntimeError, "I cannot cope with ANISOU! {0}".format(line)
#             
#             # Stop at TER
#             if line.startswith("TER"):
#                 # we write out our own TER
#                 _write_matching_residues( chainIdx, refResidues, target_residues, out )
#                 out.write("TER\n")
#                 continue
#             
#             if line.startswith("ATOM"):
#                 
#                 atom = pdb_model.PdbAtom( line )
#                 
#                 # different/first chain
#                 if atom.chainID != chain:
#                     chain = atom.chainID
#                     chainIdx+=1
#                     if chainIdx in target_residues:
#                         raise RuntimeError, "ENCOUNTERED CHAIN IN TARGET AGAIN! {0}".format( line )
#                     target_residues[ chainIdx ] = {}
#                     
#                 # We copy resSeq to make sure we don't use a reference for our index
#                 resSeq = copy.copy( atom.resSeq )
#                 
#                 # Skip any refResidues that don't match
#                 if resSeq in refResidues[ chainIdx ]:
#                 
#                     # If this is the first one add the empty tuple and reset reading
#                     if reading != resSeq:
#                         # each tuple is a list of atom objects and lines
#                         target_residues[ chainIdx ][ resSeq ] = ( [], [] )
#                         reading = resSeq
#                         
#                     target_residues[ chainIdx ][ resSeq ][0].append( atom )
#                     target_residues[ chainIdx ][ resSeq ][1].append( line )
#                     
#                 # we don't write out any atom lines as they are either not matching or 
#                 # we write out matching at the end
#                 continue
#             
#             # For time being exclude all HETATM lines
#             elif line.startswith("HETATM"):
#                 continue
#             #Endif line.startswith("ATOM")
#             
#             # Output everything else
#             out.write(line)
#             
#         # End reading loop
#         
#         # For some PDBS there is no ending TER so we need to check if we've written this out yet or not
#         if target_residues.has_key( chainIdx ):
#             _write_matching_residues( chainIdx, refResidues, target_residues, out )
#             out.write("TER\n\n")
#         
#         t.close()
#         out.close()
#         
#         return
    
    def get_info(self, inpath):
        """Read a PDB and extract as much information as possible into a PdbInfo object
        """
        
        info = pdb_model.PdbInfo()
        currentModel = None
        currentChain = -1
        
        # Go through refpdb and find which ref_residues are present
        f = open(inpath, 'r')
        line = f.readline()
        while line:
            
            # First line of title
            if line.startswith('TITLE') and not info.title:
                info.title = line[10:-1].strip()
            
            if line.startswith("REMARK"):
                
                # Resolution
                if int(line[9]) == 2:
                    line = f.readline()
                    if line.find("RESOLUTION") != -1:
                        info.resolution = float( line[25:30] )
                
                # Get solvent content                
                if int(line[7:10]) == 280:
                    
                    maxread = 5
                    # Clunky - read up to maxread lines to see if we can get the information we're after
                    # We assume the floats are at the end of the lines
                    for _ in range( maxread ):
                        line = f.readline()
                        if line.find("SOLVENT CONTENT") != -1:
                            info.solventContent = float( line.split()[-1] )
                        if line.find("MATTHEWS COEFFICIENT") != -1:
                            info.matthewsCoefficient = float( line.split()[-1] )
            #End REMARK


            if line.startswith("MODEL"):
                if currentModel:
                    # Need to make sure that we have an id if only 1 chain and none given
                    if len( currentModel.chains ) <= 1:
                        if currentModel.chains[0] == None:
                            currentModel.chains[0] = 'A'
                            
                    info.models.append( currentModel )
                    currentChain = -1
                    
                # New/first model
                currentModel = pdb_model.PdbModel()
                # Get serial
                currentModel.serial = int(line.split()[1])
            
            # Check for the first model
            if not currentModel:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    
                    # This must be the first model and there should only be one
                    currentModel = pdb_model.PdbModel()
            
            # Count chains (could also check against the COMPND line if present?)
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if line.startswith('ATOM'):
                    atom = pdb_model.PdbAtom(line)
                elif line.startswith('HETATM'):
                    atom = pdb_model.PdbHetatm(line)
            
                if atom.chainID != currentChain:    
                    # Need to check if we already have this chain for this model as a changing chain could be a sign
                    # of solvent molecules
                    if atom.chainID not in currentModel.chains:
                        currentModel.chains.append( atom.chainID )
                    currentChain = atom.chainID
            
            # Can ignore TER and ENDMDL for time being as we'll pick up changing chains anyway,
            # and new models get picked up by the models line

            line = f.readline()
            # End while loop
        
        # End of reading loop so add the last model to the list
        info.models.append( currentModel )
                    
        f.close()
        
        return info
        
    def reliable_sidechains(self, inpath=None, outpath=None ):
        """Only output non-backbone atoms for residues in the res_names list.
        """
        
        # Remove sidechains that are in res_names where the atom name is not in atom_names
        res_names = [ 'MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']
        atom_names = [ 'N', 'CA', 'C', 'O', 'CB' ]

        #   print 'Found ',each_file
        pdb_in = open( inpath, "r" )
        pdb_out = open( outpath, "w" )
        
        for pdbline in pdb_in:
            pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
            pdb_result = pdb_pattern.match(pdbline)
            
            # Check ATOM line and for residues in res_name, skip any that are not in atom names
            if pdb_result:
                pdb_result2 = re.split(pdb_pattern, pdbline)
                if pdb_result2[3] in res_names and not pdb_result2[2] in atom_names:
                    continue
            
            # Write out everything else
            pdb_out.write(pdbline)
        
        #End for
        pdb_out.close()
        pdb_in.close()
        
        return
    
    def select_residues(self, inpath=None, outpath=None, residues=None ):
        """Create a new pdb by selecting only the numbered residues from the list.
        
        Args:
        infile: path to input pdb
        outfile: path to output pdb
        residues: list of integers of the residues to keep
        
        Return:
        path to new pdb or None
        """
    
        assert inpath, outpath
        assert type(residues) == list
    
        pdb_in = open(inpath, "r")
        pdb_out = open(outpath , "w")
        
        # Loop through PDB files and create new ones that only contain the residues specified in the list
        for pdbline in pdb_in:
            pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
            pdb_result = pdb_pattern.match(pdbline)
            if pdb_result:
                pdb_result2 = re.split(pdb_pattern, pdbline )
                for i in residues : #convert to ints to comparex
        
                    if int(pdb_result2[5]) == int(i):
                        pdb_out.write(pdbline)
        
        pdb_out.close()
        
        return

    def standardise( self, inpdb, outpdb ):
        """Rename any non-standard AA, remove solvent and only keep most probably conformation.
        """
    
        tmp1 = ample_util.tmpFileName() + ".pdb" # pdbcur insists names have a .pdb suffix
        
        # Now clean up with pdbcur
        logfile = tmp1+".log"
        cmd="pdbcur xyzin {0} xyzout {1}".format( inpdb, tmp1 ).split()
        stdin="""delsolvent
    noanisou
    mostprob
    """
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
        if retcode == 0:
            # remove temporary files
            os.unlink(logfile)
        else:
            raise RuntimeError,"Error standardising pdb!"
        
        # Standardise AA names
        tmp2 = ample_util.tmpFileName() + ".pdb" # pdbcur insists names have a .pdb suffix
        self.std_residues( tmp1, tmp2 )
        
        # Strip out any remaining HETATM
        self.strip_hetatm( tmp2, outpdb )
        
        os.unlink(tmp1)
        os.unlink(tmp2) 
        
        return retcode

    def std_residues( self, pdbin, pdbout ):
        """Switch any non-standard AA's to their standard names.
        We also remove any ANISOU lines.
        """
        
        modres = [] # List of modres objects
        modres_names = {} # list of names of the modified residues keyed by chainID
        gotModel=False # to make sure we only take the first model
        reading=False # If reading structure
        
        
        pdbinf = open(pdbin,'r')
        pdboutf = open(pdbout,'w')
        
        line = True # Just for the first line
        while line:
    
            # Read in the line
            line = pdbinf.readline()
                    
            # Skip any ANISOU lines
            if line.startswith("ANISOU"):
                continue
            
            # Extract all MODRES DATA
            if line.startswith("MODRES"):
                modres.append( pdb_model.PdbModres( line ) )
                
            # Only extract the first model
            if line.startswith("MODEL"):
                if gotModel:
                    raise RuntimeError,"Found additional model! {0}".format( line )
                else:
                    gotModel=True
            
            # First time we hit coordinates we set up our data structures
            if not reading and ( line.startswith("HETATM") or line.startswith("ATOM") ):
                # There is a clever way to do this with list comprehensions but this is not it...
                for m in modres:
                    chainID = copy.copy( m.chainID )
                    if not modres_names.has_key( chainID ):
                        modres_names[ chainID ] = []
                    if m.resName not in modres_names[ chainID ]:
                        modres_names[ chainID ].append( m.resName )
                        
                # Now we're reading
                reading=True
                    
            # Switch any residue names
            if len( modres):
                if line.startswith("HETATM"):
                    
                    hetatm = pdb_model.PdbHetatm( line )
                    
                    # See if this HETATM is in the chain we are reading and one of the residues to change
                    if hetatm.resName in modres_names[ hetatm.chainID ]:
                        for m in modres:
                            if hetatm.resName == m.resName and hetatm.chainID == m.chainID:
                                # Change this HETATM to an ATOM
                                atom = pdb_model.PdbAtom().fromHetatm( hetatm )
                                # Switch residue name
                                atom.resName = m.stdRes
                                # Convert to a line
                                line = atom.toLine()+"\n"
                                break
            
            # Any HETATM have been dealt with so just process as usual
            if line.startswith("ATOM"):
                atom = pdb_model.PdbAtom( line )
                
                if atom.resName not in three2one:
                    raise RuntimeError, "Unrecognised residue! {0}".format(line)
                    
            # Output everything else
            pdboutf.write( line )
    
            # END reading loop
            
        return
    
    def strip_hetatm( self, inpath, outpath):
        """Remove all hetatoms from pdbfile"""
        o = open( outpath, 'w' )
        
        hremoved=-1
        for i, line in enumerate( open(inpath) ):
            
            # Remove EOL
            line = line.rstrip( "\n" )
            
            # Remove any HETATOM lines and following ANISOU lines
            if line.startswith("HETATM"):
                hremoved = i
                continue
            
            if line.startswith("ANISOU") and i == hremoved+1:
                continue
            
            o.write( line + "\n" )
            
        o.close()
        
        return
 
    def to_single_chain( self, inpath, outpath):
        """Condense a single-model multi-chain pdb to a single-chain pdb"""
        
        o = open( outpath, 'w' )
        
        firstChainID = None
        currentResSeq = 1 # current residue we are reading - assume it always starts from 1
        globalResSeq = 1
        for line in open(inpath):
            
            # Remove any HETATOM lines and following ANISOU lines
            if line.startswith("HETATM") or line.startswith("MODEL") or line.startswith("ANISOU"):
                raise RuntimeError,"Cant cope with the line: {0}".format( line )
            
            if line.startswith("ATOM"):
                
                changed=False
                
                atom = pdb_model.PdbAtom( line )
                
                # First atom/residue
                if not firstChainID:
                    firstChainID = atom.chainID
                
                # Change residue numbering and chainID
                if atom.chainID != firstChainID:
                    atom.chainID = firstChainID
                    changed=True
                
                # Catch each change in residue
                if atom.resSeq != currentResSeq:
                    # Change of residue
                    currentResSeq = atom.resSeq
                    globalResSeq += 1
                
                # Only change if don't match global
                if atom.resSeq != globalResSeq:
                    atom.resSeq = globalResSeq
                    changed=True
                    
                if changed:
                    line = atom.toLine()+"\n"
            
            o.write( line )
            
        o.close()
        
        return



if False:       
    #
    # Command-line handling
    #
    import argparse
    parser = argparse.ArgumentParser(description='Manipulate PDB files', prefix_chars="-")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-one_std_chain', action='store_true',
                       help='Take pdb to one model/chain that contains only standard amino acids')
    
    group.add_argument('-keep_matching', action='store_true',
                       help='keep matching atoms')
    
    parser.add_argument('-ref_file', type=str,
                       help='The reference file')
    
    parser.add_argument('input_file',
                       help='The input file - will not be altered')
    
    parser.add_argument('output_file',
                       help='The output file - will be created')
    
    if "__name__" == "__main__":
        args = parser.parse_args()
        
        # Get full paths to all files
        args.input_file = os.path.abspath( args.input_file )
        if not os.path.isfile(args.input_file):
            raise RuntimeError, "Cannot find input file: {0}".format( args.input_file )
        args.output_file = os.path.abspath( args.output_file )
        if args.ref_file:
            args.ref_file = os.path.abspath( args.ref_file )
            if not os.path.isfile(args.ref_file):
                raise RuntimeError, "Cannot find ref file: {0}".format( args.ref_file )
        
    #     if args.one_std_chain:
    #         to_1_std_chain( args.input_file, args.output_file )
    #     elif args.keep_matching:
    #         keep_matching( args.ref_file, args.input_file, args.output_file )
