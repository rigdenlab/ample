REM Run ample on windows

REM Need to add path to shelxe
set PATH=C:\Users\jmht42\Shelx\shelx64;%PATH%

%CCP4%\bin\ample.bat ^
-fasta input/toxd_.fasta ^
-mtz input/1dtx.mtz ^
-percent 50 ^
-quick_mode True ^
-nproc 6 ^
-models ..\..\tests\testfiles\models






