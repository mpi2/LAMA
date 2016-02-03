@echo off
set /p configpath=Enter (drag and drop file should work) the path to the *pheno_detect.yaml config file 
echo.

set /p projectdir=Enter (drag and drop file should) the path to the project directory (containing mutant 'inputs' folder)
echo.

python %~dp0..\pheno_detect.py -c %configpath% -p %projectdir%

cmd /k