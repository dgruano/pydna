echo "build.bat running"
"%PYTHON%" setup.py build
"%PYTHON%" setup.py install
if errorlevel 1 exit 1
