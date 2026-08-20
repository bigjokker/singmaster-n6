@echo on
REM Walk the 369 fat-image triples past prime 3.
REM Pre-registered test of Claude's size law. Not a fishing trip.
REM Predictions: s4~44 (32-58), s5~6, s6~0.8, s7~0.1, max_run 5 or 6.
REM Image clause only. Stops at kill, cell end, digit-0, or run 12.
REM Digit-0 after leaving the cell is NOT image-run length.
REM One core. Checkpoint jsonl after each k; rerun this bat to resume.
REM Minutes, not overnight. Not until-kill. Not Band II. Not E.
REM Do NOT start a second copy.
REM Requires results\fat_image_hunt.json.

cd /d "%~dp0.."
if not exist results mkdir results

if not exist results\fat_image_hunt.json (
  echo results\fat_image_hunt.json missing. Cannot walk.
  goto :fail
)

if exist results\walk_369.json (
  echo results\walk_369.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === walk 369  past prime 3  size-law test ===
echo     Source:     results\fat_image_hunt.json
echo     Checkpoint: results\walk_369.jsonl
echo     Final:      results\walk_369.json
echo     Predict:    s4~44  s5~6  s6~0.8  s7~0.1  max_run 5-6
echo.

python scripts\walk_369.py
if errorlevel 1 goto :fail

echo.
echo === walk 369 finished ===
pause
exit /b 0

:fail
echo Walk 369 was not started, or JSON already exists, or an error.
pause
exit /b 1
