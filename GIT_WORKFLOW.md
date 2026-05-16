# Git Workflow for the Local BLAST+ Prototype

Use this folder, `blast_flask_app`, as the Git repository root. Do not initialise
Git in the parent dissertation folder, because it contains large videos,
documents, downloaded BLAST+ binaries, and possible datasets.

## One-time setup

Install Git for Windows or GitHub Desktop first if `git --version` does not work.

```powershell
cd 'C:\Projects\blast_flask_app'
git init
git branch -M main
git status
```

## Suggested milestone history

Because the prototype already exists, the cleanest honest history is to commit
the current state as a set of logical milestones from this point onward. If Git
was not tracking the files when earlier changes were made, avoid pretending each
old change happened separately unless you deliberately stage file groups to
reconstruct that history.

### Option A: simplest current-state history

```powershell
git add .
git commit -m "Initial local BLAST Flask prototype"
```

Then future work can be committed in smaller pieces, such as:

```powershell
git add blast_runner.py requirements.txt templates/index.html templates/results.html smoke_test.py
git commit -m "Add Biopython SearchIO parsing"

git add blast_runner.py templates/index.html app.py
git commit -m "Add BLAST output format selection"
```

### Option B: reconstruct logical milestones from current files

If you want separate commits for the milestones already reached, stage only the
files relevant to each milestone. This is more work but gives a clearer project
story.

```powershell
git add README.md requirements.txt config.py app.py blast_runner.py smoke_test.py templates/index.html templates/results.html sample_data/
git commit -m "Create minimal Flask BLASTN prototype"

git add requirements.txt blast_runner.py smoke_test.py
git commit -m "Parse BLAST tabular and XML results with Biopython SearchIO"

git add app.py templates/index.html templates/results.html
git commit -m "Expose output format selection in the web interface"

git add .gitignore GIT_WORKFLOW.md
git commit -m "Add Git workflow and ignore generated files"
```

Only use Option B if the staged content genuinely matches the intended commit.
Check with `git diff --staged` before each commit.

## Connect to GitHub

Create an empty GitHub repository first. Do not initialise it with a README,
license, or `.gitignore`, because those files already exist locally.

Then connect this local repo to GitHub:

```powershell
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
git push -u origin main
```

If you prefer SSH:

```powershell
git remote add origin git@github.com:YOUR_USERNAME/YOUR_REPO_NAME.git
git push -u origin main
```

## Suggested GitHub issues

- Add robust FASTA validation.
- Add BLAST program selector.
- Capture runtime for each BLAST run.
- Add database registry and database creation with `makeblastdb`.
- Display command/stdout in a reproducibility panel.
- Add correctness tests against direct command-line BLAST+.
- Improve results table for clinician-friendly interpretation.
