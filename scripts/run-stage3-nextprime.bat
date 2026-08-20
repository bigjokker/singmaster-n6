@echo on
REM Stage 3: i=8 next-prime sweep k=100001..1000000.
REM Walk primes p>k until a kill. No 200/800 prime failure bucket.
REM Safety only: stop a single k if p-k exceeds 20000 (hang guard, not science).
REM One core. Expect about 8-20 hours.
REM Checkpoint jsonl: sleep/crash, run this bat again to resume.
REM Final JSON only at the end. Do not start a second copy.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\nextprime_i8_k100001-1000000.json (
  echo results\nextprime_i8_k100001-1000000.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === Stage 3  i=8  k=100001..1000000  until-killed  max-gap=20000 ===
echo     One core. Does not build m. Records q(k), q(k)-k, r(k).
echo     Checkpoint: results\nextprime_i8_k100001-1000000.jsonl
echo.

python scripts\nextprime_sweep.py --i 8 --kmin 100001 --kmax 1000000 --nprimes 0 --max-gap 20000 --json_out results\nextprime_i8_k100001-1000000.json --checkpoint results\nextprime_i8_k100001-1000000.jsonl
if errorlevel 1 goto :fail

echo.
echo === Stage 3 finished ===
pause
exit /b 0

:fail
echo Stage 3 was not started, or some k hit the 20000 gap cap / an error.
pause
exit /b 1
