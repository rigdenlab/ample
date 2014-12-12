# Run ample on windows

# Need to add path to shelxe
set PATH=C:\Users\jmht42\Shelx\shelx64;%PATH%

# Run ample
C:\Users\jmht42\Downloads\ccp4-win32-2014-11-20-0015\6.4\bin\ccp4-python ^
C:\Users\jmht42\Desktop\AMPLE\ample-dev1\bin\ample.py ^
-fasta toxd_.fasta ^
-mtz 1dtx.mtz ^
-percent 50 ^
-quick_mode True ^
-nproc 6 ^
-models_dir models



#-ensembles_dir ROSETTA_MR_0\ensembles_1




