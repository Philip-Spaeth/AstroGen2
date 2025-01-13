@echo off
setlocal

REM Pfad des Verzeichnisses, in dem die Batch-Datei liegt
set "WIN_PATH=%~dp0"

REM Entfernt den abschließenden Backslash, falls vorhanden
set "WIN_PATH=%WIN_PATH:~0,-1%"

REM Konvertiere den Windows-Pfad in einen WSL-Pfad
for /f "usebackq tokens=*" %%i in (`wsl wslpath "%WIN_PATH%"`) do set "WSL_PATH=%%i"

REM Überprüfen, ob der WSL-Pfad erfolgreich gesetzt wurde
if "%WSL_PATH%"=="" (
    echo Fehler: WSL-Pfad konnte nicht ermittelt werden.
    exit /b 1
)

REM Öffnet VSCode im WSL-Modus im Verzeichnis der Batch-Datei
wsl.exe code "%WSL_PATH%"

if errorlevel 1 (
    echo Fehler beim Öffnen von Visual Studio Code.
) else (
    echo Visual Studio Code im WSL-Modus erfolgreich gestartet.
)

endlocal