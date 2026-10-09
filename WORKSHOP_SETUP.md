# Workshop setup

Please complete these steps **before** the workshop. The slowest step (compiling
Stan) takes about 10 minutes and occasionally needs a computer restart, so
don't leave this for the morning of the session. Budget 20-30 minutes total.

If anything fails, email the organizer the output you see — don't worry about
diagnosing it yourself.

## 1. Install R and RStudio

- R (version 4.1 or newer): https://cran.r-project.org
- RStudio Desktop (free): https://posit.co/download/rstudio-desktop/

If you already have both, open RStudio and check your R version with:

```r
R.version.string
```

## 2. Install a C++ compiler

This is the step that most often goes wrong, so follow it carefully for your
operating system. You won't interact with the compiler directly — it runs in
the background to build the statistical model — but it has to be present.

**Windows:** Install Rtools, matching your R version, from
https://cran.r-project.org/bin/windows/Rtools/. Use the default options
during installation.

**Mac:** Open Terminal (Applications > Utilities > Terminal) and run:

```
xcode-select --install
```

Click "Install" in the dialog that pops up. This can take several minutes.

**Linux:** Install `g++` with your package manager, e.g. on Ubuntu/Debian:

```
sudo apt install g++
```

## 3. Install Stan (cmdstanr)

In RStudio, run:

```r
install.packages("cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

cmdstanr::install_cmdstan()   # one-time, ~10 minutes
```

## 4. Install the countsynth package

```r
install.packages("remotes")
remotes::install_github("afranks86/countsynth", dependencies = TRUE)
```

## 5. Run the setup check

Download
[`check_workshop_setup.R`](scripts/check_workshop_setup.R)
and run it in RStudio (open it and click "Source"), or from a terminal:

```
Rscript check_workshop_setup.R
```

It checks your R version, compiler, Stan installation, and the countsynth
package, then runs an end-to-end test to confirm everything works together.
**The last step takes several minutes** (5-10 is normal, and the very first
run is slower because it has to compile the model) — let it finish. You
should see:

```
All checks passed. You're ready for the workshop.
```

If you see `FAILED` next to any step, the script prints what to do next. If
you're still stuck after trying that, email the organizer the full output —
screenshots are fine.
