@echo on
REM Overnight Band I census: exhaustive image runs on fat non-zero cells.
REM NONE + PART-lower, P_hi > 1e6, g_max > 25000. Largest g first.
REM Prize window: (3,1) PART-lower, rho~0.176, g~475000.
REM First three primes only. Looking for a triple (or max run still 2).
REM One core. Checkpoint jsonl after each region; rerun this bat to resume.
REM Expect most of a night. Not until-kill. Not Band II. Not E.
REM Do NOT start a second copy.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\fat_image_hunt.json (
  echo results\fat_image_hunt.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === OVERNIGHT fat image hunt  exhaust g^>25000  P_hi^>1e6 ===
echo     Checkpoint: results\fat_image_hunt.jsonl
echo.

python scripts\fat_image_hunt.py
if errorlevel 1 goto :fail

echo.
echo === fat image hunt finished ===
pause
exit /b 0

:fail
echo Fat image hunt was not started, or JSON already exists, or an error.
pause
exit /b 1
