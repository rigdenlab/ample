REM Toxd test case

%CCP4%\bin\ample.bat ^
-fasta input\toxd_.fasta ^
-mtz input\1dtx.mtz ^
-percent 50 ^
-quick_mode True ^
-nproc 6 ^
-models ..\..\tests\testfiles\models ^
-show_gui True

