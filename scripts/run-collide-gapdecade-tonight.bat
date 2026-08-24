@echo on
REM Skipped collide decade, TONIGHT pack: l=20,19,18,17,16,9,15,8,14,13.
REM The 61/101-digit hole D5 left unscanned. Engine computes m-range;
REM does not overlap recorded collide_*.json (max-m = old start - 1).
REM ~1.9 h serial, ~15-25 min at 8 workers. 127 pairs.
REM Skip-if-json-exists: rerun this bat to resume. Do NOT start a second copy.
REM Not nearby. Not l=10 (20 h). Not (3,5)/(4,5). Not i=10.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\collide_gapdecade_tonight.lock (
  echo results\collide_gapdecade_tonight.lock exists. A copy is running, or a crash left the lock.
  echo Delete the lock only if no python collide is running, then rerun to resume.
  goto :fail
)

echo.
echo === TONIGHT collide gap-decade  l=20..13  workers=8 ===
echo     Driver:     scripts\collide_gapdecade.py --pack tonight
echo     Outputs:    results\collide_gapdecade_k*_l*.json
echo     Resume:     rerun this bat; existing json are skipped
echo.

echo running> results\collide_gapdecade_tonight.lock
python scripts\collide_gapdecade.py --pack tonight --workers 8
set ERR=%ERRORLEVEL%
del results\collide_gapdecade_tonight.lock
if not "%ERR%"=="0" goto :fail

echo.
echo === tonight pack finished ===
pause
exit /b 0

:fail
echo Tonight collide pack was not started, or a pair failed.
pause
exit /b 1
