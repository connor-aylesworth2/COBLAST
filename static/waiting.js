(function () {
  // This file gives long-running POST forms an immediate full-screen status
  // page, then lets the normal form submission continue in the background.
  const DEFAULT_MESSAGES = [
    "Preparing the local job.",
    "Working locally on this machine.",
    "Large files can take several minutes.",
    "Elapsed time is still increasing while the job is active.",
  ];

  let overlay = null;
  let elapsedTimer = null;
  let messageTimer = null;

  function formatElapsed(totalSeconds) {
    // Render elapsed time in the same HH:MM:SS style as remote BLAST pages.
    const hours = Math.floor(totalSeconds / 3600);
    const minutes = Math.floor((totalSeconds % 3600) / 60);
    const seconds = totalSeconds % 60;

    return `${String(hours).padStart(2, "0")}:${String(minutes).padStart(2, "0")}:${String(seconds).padStart(2, "0")}`;
  }

  function formatDateTime(date) {
    // Locale-aware formatting keeps the display familiar on different machines.
    return date.toLocaleString(undefined, {
      weekday: "short",
      year: "numeric",
      month: "short",
      day: "2-digit",
      hour: "2-digit",
      minute: "2-digit",
      second: "2-digit",
    });
  }

  function localRequestId() {
    // Local jobs do not have real NCBI RIDs, so generate a readable stand-in.
    const suffix = Math.random().toString(36).slice(2, 8).toUpperCase();
    return `LOCAL-${Date.now().toString(36).toUpperCase()}-${suffix}`;
  }

  function messagesForForm(form) {
    // Forms can customize progress messages with a pipe-separated data attribute.
    const rawMessages = form.dataset.waitMessages || "";
    const messages = rawMessages
      .split("|")
      .map((message) => message.trim())
      .filter(Boolean);
    return messages.length ? messages : DEFAULT_MESSAGES;
  }

  function ensureOverlay() {
    // Lazily build the overlay once, then reuse it for every waiting form.
    if (overlay) {
      return overlay;
    }

    overlay = document.createElement("div");
    // role/status and aria-live announce progress changes to assistive tech.
    overlay.className = "working-overlay";
    overlay.setAttribute("role", "status");
    overlay.setAttribute("aria-live", "polite");
    overlay.innerHTML = `
      <div class="working-page">
        <header class="working-ncbi-header">
          <div class="working-ncbi-inner">
            <div class="working-logo-text">
              <strong>COBLAST+</strong>
              <span>Clinician-Oriented BLAST+</span>
            </div>
          </div>
        </header>
        <div class="working-nav-strip">
          <div class="working-nav-inner">
            <div class="working-breadcrumb">
              <strong>COBLAST+</strong> &raquo; local job &raquo; <span class="working-rid-breadcrumb">RID-LOCAL</span>
            </div>
            <nav class="working-links" aria-label="Waiting page navigation">
              <span>Home</span>
              <span>Recent Results</span>
              <span>Saved Strategies</span>
              <span>Help</span>
            </nav>
          </div>
        </div>
        <main class="working-status-page">
          <h3 class="working-job-title">Job Title: <span>Working</span></h3>

          <table class="working-status-table">
            <tbody>
              <tr>
                <th>Request ID</th>
                <td class="working-request-id">LOCAL</td>
              </tr>
              <tr>
                <th>Status</th>
                <td class="working-status">Searching</td>
              </tr>
              <tr>
                <th>Submitted at</th>
                <td class="working-submitted-at"></td>
              </tr>
              <tr>
                <th>Current time</th>
                <td class="working-current-time"></td>
              </tr>
              <tr>
                <th>Time since submission</th>
                <td class="working-elapsed">00:00:00</td>
              </tr>
            </tbody>
          </table>

          <p class="working-message">This page will be automatically updated in <strong>2</strong> seconds until the local job is done.</p>
          <p class="working-detail">Starting local job.</p>
        </main>
      </div>
    `;
    document.body.appendChild(overlay);
    return overlay;
  }

  function showWaitingScreen(form) {
    // Populate the overlay from the submitted form's data attributes.
    const currentOverlay = ensureOverlay();
    const title = form.dataset.waitTitle || "Working";
    const messages = messagesForForm(form);
    const requestId = localRequestId();
    const submittedAt = new Date();
    const breadcrumbNode = currentOverlay.querySelector(".working-rid-breadcrumb");
    const titleNode = currentOverlay.querySelector(".working-job-title span");
    const requestIdNode = currentOverlay.querySelector(".working-request-id");
    const statusNode = currentOverlay.querySelector(".working-status");
    const submittedAtNode = currentOverlay.querySelector(".working-submitted-at");
    const currentTimeNode = currentOverlay.querySelector(".working-current-time");
    const elapsedNode = currentOverlay.querySelector(".working-elapsed");
    const detailNode = currentOverlay.querySelector(".working-detail");
    let messageIndex = 0;

    clearInterval(elapsedTimer);
    clearInterval(messageTimer);

    // Reset every changing field so repeat submissions start from a clean state.
    breadcrumbNode.textContent = `RID-${requestId}`;
    titleNode.textContent = title;
    requestIdNode.textContent = requestId;
    statusNode.textContent = "Searching";
    submittedAtNode.textContent = formatDateTime(submittedAt);
    currentTimeNode.textContent = formatDateTime(submittedAt);
    elapsedNode.textContent = "00:00:00";
    detailNode.textContent = messages[messageIndex];
    currentOverlay.classList.add("is-visible");
    document.body.classList.add("is-waiting");

    elapsedTimer = setInterval(() => {
      // Keep elapsed/current time moving while the Flask request is running.
      const now = new Date();
      const elapsedSeconds = Math.floor((now.getTime() - submittedAt.getTime()) / 1000);
      currentTimeNode.textContent = formatDateTime(now);
      elapsedNode.textContent = formatElapsed(elapsedSeconds);
    }, 1000);

    messageTimer = setInterval(() => {
      // Cycle through helpful messages so the user knows the page is alive.
      messageIndex = (messageIndex + 1) % messages.length;
      detailNode.textContent = messages[messageIndex];
    }, 9000);
  }

  function attachWaitingScreens() {
    // Any form with data-wait-title opts into the waiting overlay behavior.
    const forms = document.querySelectorAll("form[data-wait-title]");
    for (const form of forms) {
      form.addEventListener("submit", (event) => {
        if (form.dataset.waitSubmitted === "true") {
          return;
        }

        event.preventDefault();
        form.dataset.waitSubmitted = "true";
        showWaitingScreen(form);

        // Disable submit buttons to avoid duplicate BLAST jobs.
        for (const button of form.querySelectorAll("button[type='submit'], input[type='submit']")) {
          button.disabled = true;
        }

        window.setTimeout(() => {
          // Calling the prototype submit avoids recursively firing this handler.
          HTMLFormElement.prototype.submit.call(form);
        }, 80);
      });
    }
  }

  window.addEventListener("pageshow", () => {
    // Browser back/forward cache can restore pages, so always hide stale overlay.
    document.body.classList.remove("is-waiting");
    if (overlay) {
      overlay.classList.remove("is-visible");
    }
    clearInterval(elapsedTimer);
    clearInterval(messageTimer);
  });

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attachWaitingScreens);
  } else {
    attachWaitingScreens();
  }
})();
