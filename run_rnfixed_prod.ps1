# run_rnfixed_prod.ps1 -- resumable production driver (eps_rn_me fixed at 0.0015, mh_jscale=0.01,
# prior_variance jumping): 10 chunks x 5000 draws x 2 chains = 50,000 draws/chain. Runs each chunk
# in its own MATLAB process at BELOW-NORMAL priority. Progress is stored in
# rnfixed_prod_chunks_done.txt, so re-running this script after an interruption resumes at the
# next chunk. Estimated ~12h total based on the 8000-draw diagnostic's pace (118.6 min for 8000
# draws/chain), well inside the 36h+ window the computer will stay awake.
# To stop: create a file named STOP_RNFIXED_PROD in this folder (the driver stops after the
# current chunk), or end the MATLAB/powershell processes (an interrupted chunk should be reported
# before resuming).
$dir   = "C:\Users\gawater\Documents\GitHub\Shotgun-2"
$total = 10
$stateFile = Join-Path $dir "rnfixed_prod_chunks_done.txt"
Set-Location $dir
$done = 0
if (Test-Path $stateFile) { $done = [int](Get-Content $stateFile -Raw).Trim() }
while ($done -lt $total) {
    if (Test-Path (Join-Path $dir "STOP_RNFIXED_PROD")) { "stop file found; stopping after $done chunks" | Add-Content (Join-Path $dir "rnfixed_prod_driver.log"); break }
    $loadflag = 0; if ($done -gt 0) { $loadflag = 1 }
    $log = Join-Path $dir ("rnfixed_prod_chunk{0:D2}.log" -f ($done + 1))
    Remove-Item $log -ErrorAction SilentlyContinue
    "$(Get-Date -Format s) starting chunk $($done + 1) of $total (loadflag=$loadflag)" | Add-Content (Join-Path $dir "rnfixed_prod_driver.log")
    $cmd = "loadflag=$loadflag; run('C:/Users/gawater/Documents/GitHub/Shotgun-2/run_rnfixed_chunk.m')"
    $argString = "-batch `"$cmd`" -logfile `"$log`""
    $p = Start-Process -FilePath "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -ArgumentList $argString -WindowStyle Hidden -PassThru
    while (-not $p.HasExited) {
        Start-Sleep -Seconds 20
        Get-Process MATLAB -ErrorAction SilentlyContinue | Where-Object { $_.StartTime -ge $p.StartTime.AddSeconds(-5) } | ForEach-Object { try { $_.PriorityClass = 'BelowNormal' } catch {} }
    }
    $ok = (Test-Path $log) -and ((Get-Content $log -Raw) -match 'CHUNK FINISHED')
    if (-not $ok) { "$(Get-Date -Format s) chunk $($done + 1) did NOT finish cleanly; stopping. See $log" | Add-Content (Join-Path $dir "rnfixed_prod_driver.log"); break }
    $done++
    "$done" | Set-Content $stateFile
    "$(Get-Date -Format s) chunk $done finished" | Add-Content (Join-Path $dir "rnfixed_prod_driver.log")
}
"$(Get-Date -Format s) driver exiting with $done of $total chunks done" | Add-Content (Join-Path $dir "rnfixed_prod_driver.log")
