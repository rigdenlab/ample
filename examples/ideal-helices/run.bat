REM Windows run script

REM Need to add path to shelxe
set PATH=C:\Users\jmht42\Shelx\shelx64;%PATH%

%CCP4%\bin\ample.bat ^
-fasta input\2OVC.fasta ^
-mtz input\2OVC-cad.mtz ^
-use_shelxe True ^
-ideal_helices True ^
-nproc 8 ^
-show_gui True
