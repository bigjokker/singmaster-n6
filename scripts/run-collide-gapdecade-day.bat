@echo on
REM Skipped collide decade, DAY pack: l=12,7,11,6.
REM Run AFTER tonight pack succeeds. ~7 h at 8 workers.
REM Not l=10 (~20 h at 8 workers) -- that is a later overnight.
REM Skip-if-json-exists. Do NOT start a second copy. Do NOT start if tonight still running.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\collide_gapdecade_tonight.lock (
  echo tonight pack still running. Do not overlap.
  goto :fail
)
if exist results\collide_gapdecade_day.lock (
  echo results\collide_gapdecade_day.lock exists. A copy is running, or a crash left the lock.
  goto :fail
)

echo.
echo === DAY collide gap-decade  l=12,7,11,6  workers=8 ===
echo     Driver:     scripts\collide_gapdecade.py --pack day
echo     Expect:     ~7 hours at 8 workers
echo     Resume:     rerun this bat; existing json are skipped
echo.

echo running> results\collide_gapdecade_day.lock
python scripts\collide_gapdecade.py --pack day --workers 8
set ERR=%ERRORLEVEL%
del results\collide_gapdecade_day.lock
if not "%ERR%"=="0" goto :fail

echo.
echo === day pack finished ===
pause
exit /b 0

:fail
echo Day collide pack was not started, or a pair failed.
pause
exit /b 1
