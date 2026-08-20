@echo on
REM Stage 2: i=8 next-prime sweep k=10001..100000, up to 200 primes after each k.
REM One core. Expect about 5-15 minutes (worst case ~30 min if r grows).
REM JSON is written only when the sweep ends.
REM See docs/modular-spec.txt. Do not start a second copy.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\nextprime_i8_k10001-100000.json (
  echo results\nextprime_i8_k10001-100000.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === Stage 2  i=8  k=10001..100000  nprimes=200 ===
echo     One core. Next-prime obstruction only. Does not build m.
echo.

python scripts\nextprime_sweep.py --i 8 --kmin 10001 --kmax 100000 --nprimes 200 --json_out results\nextprime_i8_k10001-100000.json
if errorlevel 1 goto :fail

echo.
echo === Stage 2 finished ===
pause
exit /b 0

:fail
echo Stage 2 was not started, or the sweep reported unkilled columns / an error.
pause
exit /b 1
