@echo on
REM Harvest i=9 witnesses in this playground (bisect first_live_after ~72s).
REM The live harvest would take ~2.3h on the old linear scan.
REM
REM Run AFTER Desktop\Singmaster\results\i9_sweep.json exists.
REM Do NOT run while the live CMD is still in python -- the jsonl is
REM still being appended and a partial build is not a certificate.
REM Copies the jsonl here, does not delete the live file.
REM
REM 2026-08-23: THIS ALREADY RAN. results\i9_witness.npz in this clone is
REM complete -- 35,522,326 rows, n_unresolved=0, k=2..80 and the four Lucas
REM columns filled. A full clean re-run would rebuild and re-fill to the same
REM table, so the risk is not the happy path: it is an ABORTED or failing
REM re-run replacing a complete table with the unfilled intermediate the
REM build step writes, since build and fill are two separate commands here.
REM The bat therefore refuses when the table is already present, in the same
REM spirit as every other bat refusing when its json exists.
REM Set HARVEST_FORCE=1 to override.

set LIVE=C:\Users\wwwsa\Desktop\Singmaster
cd /d "%~dp0.."

if exist results\i9_witness.npz (
  if not "%HARVEST_FORCE%"=="1" (
    echo results\i9_witness.npz already exists -- i=9 is already harvested here.
    echo Nothing to do. An aborted re-run would leave the unfilled build output.
    echo Set HARVEST_FORCE=1 only if you really mean to rebuild it.
    goto :fail
  )
)

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
REM The ledger reports TWO facts and exits 1 unless BOTH hold for every
REM member: coverage complete, and a sweep certificate naming that table.
REM i=8 and i=9 are legitimately UNBOUND -- their runs ended clean=false, so
REM their records carry certificate=null, which is the honest state after the
REM k=1021 repair. Exit 1 for that reason alone is NOT a harvest failure, so
REM do not test errorlevel here. Fail on a real coverage gap (GAPS PRESENT)
REM or on a member that could not be read (ERROR) -- both are printed.
python scripts\coverage_ledger.py --json_out results\coverage_ledger.json > "%TEMP%\i9_ledger.out" 2>&1
type "%TEMP%\i9_ledger.out"
findstr /C:"GAPS PRESENT" "%TEMP%\i9_ledger.out" >nul && goto :ledger_gap
findstr /C:"ERROR" "%TEMP%\i9_ledger.out" >nul && goto :ledger_gap
findstr /C:"columns witnessed" "%TEMP%\i9_ledger.out" >nul || goto :ledger_gap

echo.
echo === verify sample (witnesses are actually kills; 5000 of ~28M) ===
python scripts\witness.py verify --file results\i9_witness.npz --sample 5000
if errorlevel 1 goto :fail

echo.
echo === i=9 playground harvest done ===
echo     KEEP the live jsonl. coverage is complete; 5000 certificates checked.
echo     full verify of 28M is a later job, not this bat.
echo     i=9 reads "coverage complete, UNBOUND" in the ledger: complete
echo     coverage, no certificate naming the table. That is expected here.
pause
exit /b 0

:ledger_gap
echo.
echo the coverage ledger reported a real gap (missing/extra columns, or a
echo member it could not read). That is not the UNBOUND state -- read the
echo output above before doing anything else.
pause
exit /b 1

:fail
echo harvest was not run, or a step failed.
pause
exit /b 1
