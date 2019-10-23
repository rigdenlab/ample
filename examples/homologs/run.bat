REM Windows run script

REM Need to add path to shelxe
set PATH=C:\Users\jmht42\Shelx\shelx64;%PATH%

%CCP4%\bin\ample.bat \
-fasta input\3DCY.fasta ^
-mtz input\3dcy-sf.mtz ^
-homologs True ^
-models .  ^
-nproc 8 ^
-show_gui True

# Add below for running in benchmark mode
#-mustang_exe C:\opt\MUSTANG_v3.2.2\bin\mustang-3.2.1 ^
#-alignment_file input\testthree.afasta ^
#-native_pdb input\1DCY.pdb ^
