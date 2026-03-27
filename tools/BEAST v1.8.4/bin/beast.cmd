@echo off
setlocal
if ""=="%BEAST%" set BEAST=%~dp0%..
set BEAST_LIB=%BEAST%\lib
java -Xms64m -Xmx256m -Djava.library.path="%BEAST_LIB%" -Dbeast.plugins.dir="%BEAST_LIB%\plugins" -jar "%BEAST_LIB%/beast.jar" %*
