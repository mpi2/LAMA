set /p configpath=Please enter the path to the lama config file

python %~dp0..\lama.py -c %configpath%

cmd /k