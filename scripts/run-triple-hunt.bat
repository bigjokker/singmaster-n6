@echo on
REM NONE-window triple hunt on two-digit cells with P_hi > 1e6.
REM Does NOT redo the 1e5-1e6 census (already max run 2, 0 triples).
REM Exhaustive under-Z k only if g_max <= 25000.
REM Fat slabs (4,1)/(7,2) get three canonical k only — not O(g^2).
REM Expect tens of minutes, not a second Stage 3.
REM Run AFTER stragglers, not instead. Do not start a second copy.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\triple_hunt_p1e6-K.json (
  echo results\triple_hunt_p1e6-K.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === triple hunt  NONE cells P_hi^>1e6  exhaust g^<=25000 ===
echo.

python scripts\triple_hunt.py
if errorlevel 1 goto :fail

echo.
echo === triple hunt finished ===
pause
exit /b 0

:fail
echo Triple hunt was not started, or JSON already exists, or an error.
pause
exit /b 1
