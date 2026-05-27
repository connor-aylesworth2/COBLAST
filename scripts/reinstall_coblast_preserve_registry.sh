#!/usr/bin/env bash
set -Eeuo pipefail

# Reinstall COBLAST from GitHub while preserving the existing SQLite registry.
#
# Typical use on the Linux server:
#   bash scripts/reinstall_coblast_preserve_registry.sh \
#     --install-dir "$HOME/COBLAST_1.00/COBLAST-"
#
# The script intentionally requires --install-dir so it never guesses which
# folder should be replaced.

REMOTE_URL="https://github.com/connor-aylesworth2/COBLAST-.git"
BRANCH="main"
INSTALL_DIR=""
REGISTRY_PATH=""
BACKUP_ROOT="${COBLAST_REINSTALL_LOG_DIR:-$HOME/COBLAST_reinstall_logs}"
DELETE_OLD=true
INSTALL_REQUIREMENTS=false

usage() {
  cat <<'USAGE'
Usage:
  reinstall_coblast_preserve_registry.sh --install-dir PATH [options]

Required:
  --install-dir PATH        Existing COBLAST checkout to replace.

Options:
  --registry-path PATH      Registry to preserve. Defaults to:
                            PATH/instance/database_registry.sqlite
  --backup-root PATH        Folder for reinstall logs/backups.
                            Defaults to ~/COBLAST_reinstall_logs
  --remote-url URL          Git remote to clone.
                            Defaults to connor-aylesworth2/COBLAST-
  --branch NAME             Branch to clone. Defaults to main.
  --keep-old                Keep the archived old checkout after success.
  --install-requirements    Run python3 -m pip install -r requirements.txt
                            after cloning. Activate conda first if using this.
  -h, --help                Show this help text.

What it does:
  1. Moves the old database_registry.sqlite into a timestamped backup folder.
  2. Moves the old COBLAST checkout aside.
  3. Clones the latest COBLAST code from GitHub into the original path.
  4. Restores the preserved SQLite registry into the new checkout.
  5. Deletes the archived old checkout after success, unless --keep-old is used.
USAGE
}

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

require_arg() {
  local option="$1"
  local value="${2:-}"
  [[ -n "$value" && "$value" != --* ]] || die "$option requires a value."
}

resolve_path() {
  python3 - "$1" <<'PY'
from pathlib import Path
import sys

print(Path(sys.argv[1]).expanduser().resolve())
PY
}

ensure_under_home() {
  local path="$1"
  case "$path" in
    "$HOME_DIR"/*) ;;
    *) die "Refusing to modify a path outside your home directory: $path" ;;
  esac
}

ensure_safe_tree_target() {
  local path="$1"
  ensure_under_home "$path"
  [[ "$path" != "$HOME_DIR" ]] || die "Refusing to modify your home directory."
  [[ "$path" != "/" ]] || die "Refusing to modify /."
}

remove_old_archive() {
  local path="$1"
  ensure_safe_tree_target "$path"
  [[ -f "$path/run_COBLAST.py" && -f "$path/app.py" ]] || {
    die "Refusing to remove $path because it does not look like a COBLAST checkout."
  }
  rm -rf -- "$path"
}

rollback_needed=false
registry_moved=false
old_moved=false
BACKUP_DIR=""
REGISTRY_BACKUP=""
OLD_ARCHIVE=""

rollback() {
  local exit_code=$?
  if [[ "$rollback_needed" != true ]]; then
    exit "$exit_code"
  fi

  log "Reinstall failed. Attempting to roll back to the previous checkout."

  if [[ "$old_moved" == true && -d "$OLD_ARCHIVE" ]]; then
    if [[ -d "$INSTALL_DIR" ]]; then
      ensure_safe_tree_target "$INSTALL_DIR"
      rm -rf -- "$INSTALL_DIR"
    fi
    mkdir -p "$OLD_ARCHIVE/instance"
    if [[ -f "$REGISTRY_BACKUP" ]]; then
      cp -p "$REGISTRY_BACKUP" "$OLD_ARCHIVE/instance/database_registry.sqlite"
    fi
    mv "$OLD_ARCHIVE" "$INSTALL_DIR"
    log "Rollback restored the old checkout at: $INSTALL_DIR"
  elif [[ "$registry_moved" == true && -f "$REGISTRY_BACKUP" ]]; then
    mkdir -p "$(dirname "$REGISTRY_PATH")"
    cp -p "$REGISTRY_BACKUP" "$REGISTRY_PATH"
    log "Rollback restored the registry at: $REGISTRY_PATH"
  fi

  log "The registry backup remains at: $REGISTRY_BACKUP"
  exit "$exit_code"
}

trap rollback ERR

while [[ $# -gt 0 ]]; do
  case "$1" in
    --install-dir)
      require_arg "$1" "${2:-}"
      INSTALL_DIR="$2"
      shift 2
      ;;
    --registry-path)
      require_arg "$1" "${2:-}"
      REGISTRY_PATH="$2"
      shift 2
      ;;
    --backup-root)
      require_arg "$1" "${2:-}"
      BACKUP_ROOT="$2"
      shift 2
      ;;
    --remote-url)
      require_arg "$1" "${2:-}"
      REMOTE_URL="$2"
      shift 2
      ;;
    --branch)
      require_arg "$1" "${2:-}"
      BRANCH="$2"
      shift 2
      ;;
    --keep-old)
      DELETE_OLD=false
      shift
      ;;
    --install-requirements)
      INSTALL_REQUIREMENTS=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

[[ -n "$INSTALL_DIR" ]] || {
  usage
  die "Missing required --install-dir PATH."
}

command -v python3 >/dev/null 2>&1 || die "python3 is required for path handling."
command -v git >/dev/null 2>&1 || die "git is required to clone COBLAST."

HOME_DIR="$(resolve_path "$HOME")"
INSTALL_DIR="$(resolve_path "$INSTALL_DIR")"
BACKUP_ROOT="$(resolve_path "$BACKUP_ROOT")"

if [[ -z "$REGISTRY_PATH" ]]; then
  REGISTRY_PATH="$INSTALL_DIR/instance/database_registry.sqlite"
else
  REGISTRY_PATH="$(resolve_path "$REGISTRY_PATH")"
fi

ensure_safe_tree_target "$INSTALL_DIR"
ensure_under_home "$BACKUP_ROOT"

[[ -d "$INSTALL_DIR" ]] || die "Install directory does not exist: $INSTALL_DIR"
[[ -f "$INSTALL_DIR/run_COBLAST.py" && -f "$INSTALL_DIR/app.py" ]] || {
  die "Install directory does not look like a COBLAST checkout: $INSTALL_DIR"
}
[[ -f "$REGISTRY_PATH" ]] || die "Registry file does not exist: $REGISTRY_PATH"

TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
BACKUP_DIR="$BACKUP_ROOT/reinstall_$TIMESTAMP"
REGISTRY_BACKUP="$BACKUP_DIR/database_registry.sqlite"
OLD_ARCHIVE="${INSTALL_DIR}.old.$TIMESTAMP"

mkdir -p "$BACKUP_DIR"
exec > >(tee -a "$BACKUP_DIR/reinstall.log") 2>&1
rollback_needed=true

log "Starting COBLAST reinstall."
log "Install directory: $INSTALL_DIR"
log "Registry source:   $REGISTRY_PATH"
log "Backup folder:     $BACKUP_DIR"
log "Remote URL:        $REMOTE_URL"
log "Branch:            $BRANCH"

[[ ! -e "$OLD_ARCHIVE" ]] || die "Archive path already exists: $OLD_ARCHIVE"

log "Moving registry into backup folder."
mv "$REGISTRY_PATH" "$REGISTRY_BACKUP"
registry_moved=true

log "Moving old checkout aside."
mv "$INSTALL_DIR" "$OLD_ARCHIVE"
old_moved=true

log "Cloning latest COBLAST checkout."
git clone --branch "$BRANCH" --single-branch "$REMOTE_URL" "$INSTALL_DIR"

log "Restoring preserved registry into the new checkout."
mkdir -p "$INSTALL_DIR/instance"
cp -p "$REGISTRY_BACKUP" "$INSTALL_DIR/instance/database_registry.sqlite"

if [[ "$INSTALL_REQUIREMENTS" == true ]]; then
  log "Installing Python requirements with python3 -m pip."
  (cd "$INSTALL_DIR" && python3 -m pip install -r requirements.txt)
fi

NEW_COMMIT="$(git -C "$INSTALL_DIR" rev-parse --short HEAD)"

if [[ "$DELETE_OLD" == true ]]; then
  log "Deleting archived old checkout."
  remove_old_archive "$OLD_ARCHIVE"
else
  log "Keeping archived old checkout at: $OLD_ARCHIVE"
fi

rollback_needed=false

log "Reinstall complete."
log "Registry backup kept at: $REGISTRY_BACKUP"
log "New checkout commit: $NEW_COMMIT"
log "Next: activate your COBLAST environment, export BLAST_BIN if needed, and start the app."
