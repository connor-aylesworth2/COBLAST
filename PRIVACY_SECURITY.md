# Privacy and Security Notes

This prototype is designed for local testing of a clinician-friendly BLAST+
interface. It is not a deployed clinical service and it is not validated for
diagnostic use.

## Local-only behavior

- The Flask development server is configured to bind to `127.0.0.1`.
- Remote BLAST is disabled.
- The backend rejects generated BLAST commands containing `-remote`.
- Query FASTA files, result files, and the SQLite database registry are stored
  locally under the application workspace, primarily in `instance\`.

## Tester responsibilities

- Use toy, public, synthetic, or otherwise approved data during external
  testing.
- Do not upload patient-identifiable, confidential, controlled-access, or
  unpublished sequence data to GitHub issues, pull requests, chat tools, or
  screenshots.
- Do not expose the Flask development server to a network without a separate
  security review.
- Do not rely on prototype output for clinical diagnosis, treatment decisions,
  or patient management.

## Current security boundary

This project wraps local BLAST+ executables through a constrained Python backend.
It is intended to reduce accidental command-line misuse, not to provide a
hardened multi-user web service. Before any networked or clinical deployment,
the project would need authentication, authorization, audit logging, stronger
input limits, server hardening, and formal data governance review.
