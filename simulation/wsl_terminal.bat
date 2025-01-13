@echo off
REM Hole den Pfad des aktuellen Verzeichnisses
set current_dir=%~dp0

REM Wandle den Pfad in das WSL-Format um
for /f %%i in ('wsl wslpath "%current_dir%"') do set wsl_dir=%%i

REM Öffne das WSL-Terminal in diesem Verzeichnis
wsl -d Ubuntu -e bash -c "cd '%wsl_dir%' && exec bash"