echo off
set LOCALHOST=%COMPUTERNAME%
set REMOTE_SH="C:\PROGRA~1\ANSYSI~1\v171\fluent/ntbin/win64/rsh.exe"
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v171\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v171\fluent\ntbin\win64\tell.exe" Enmodes1 60284 CLEANUP_EXITING
if /i "%LOCALHOST%"=="Enmodes1" (%KILL_CMD% 9696) else (%REMOTE_SH% Enmodes1 %KILL_CMD% 9696)
del "D:\Sonstiges\OwnTutorials\Multiscale\Fluent\cleanup-fluent-Enmodes1-4128.bat"
