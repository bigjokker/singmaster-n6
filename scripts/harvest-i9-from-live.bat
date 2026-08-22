@echo on
REM Harvest i=9 witnesses in this playground (bisect first_live_after ~72s).
REM The live harvest would take ~2.3h on the old linear scan.
REM
REM Run AFTER Desktop\Singmaster\results\i9_sweep.json exists.
REM Do NOT run while the live CMD is still in python -- the jsonl is
REM still being appended and a partial build is not a certificate.
REM Copies the jsonl here, does not delete the live file.

set LIVE=C:\Users\wwwsa\Desktop\Singmaster
cd /d "%~dp0.."

if not exist "%LIVE%\results\i9_sweep.json" (
  echo live i9_sweep.json is missing. i=9 is not finished.
  echo wait until the live CMD prints "i=9 sweep finished".
  goto :fail
)
if not exist "%LIVE%\results\i9_sweep.jsonl" (
  echo live i9_sweep.jsonl is missing. Cannot build witnesses.
  goto :fail
)

if not exist results mkdir results
echo.
echo === copy live i=9 artifacts into the playground ===
copy /Y "%LIVE%\results\i9_sweep.jsonl" results\i9_sweep.jsonl
if errorlevel 1 goto :fail
copy /Y "%LIVE%\results\i9_sweep.json" results\i9_sweep.json
if errorlevel 1 goto :fail

echo.
echo === build i=9 witness table (playground kernel) ===
python scripts\witness.py build --i 9 --checkpoint results\i9_sweep.jsonl --out results\i9_witness.npz
if errorlevel 1 goto :fail

echo.
echo === fill k=2..80 from the modular engine ===
python scripts\witness.py fill --i 9 --file results\i9_witness.npz
if errorlevel 1 goto :fail

echo.
echo === coverage ledger (every claimed column has a witness) ===
python scripts\coverage_ledger.py --json_out results\coverage_ledger.json
if errorlevel 1 goto :fail

echo.
echo === verify sample (witnesses are actually kills; 5000 of ~28M) ===
python scripts\witness.py verify --file results\i9_witness.npz --sample 5000
if errorlevel 1 goto :fail

echo.
echo === i=9 playground harvest done ===
echo     KEEP the live jsonl. coverage is complete; 5000 certificates checked.
echo     full verify of 28M is a later job, not this bat.
pause
exit /b 0

:fail
echo harvest was not run, or a step failed.
pause
exit /b 1
