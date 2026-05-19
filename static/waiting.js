(function () {
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
    const hours = Math.floor(totalSeconds / 3600);
    const minutes = Math.floor((totalSeconds % 3600) / 60);
    const seconds = totalSeconds % 60;

    if (hours > 0) {
      return `${hours}:${String(minutes).padStart(2, "0")}:${String(seconds).padStart(2, "0")}`;
    }
    return `${minutes}:${String(seconds).padStart(2, "0")}`;
  }

  function messagesForForm(form) {
    const rawMessages = form.dataset.waitMessages || "";
    const messages = rawMessages
      .split("|")
      .map((message) => message.trim())
      .filter(Boolean);
    return messages.length ? messages : DEFAULT_MESSAGES;
  }

  function ensureOverlay() {
    if (overlay) {
      return overlay;
    }

    overlay = document.createElement("div");
    overlay.className = "working-overlay";
    overlay.setAttribute("role", "status");
    overlay.setAttribute("aria-live", "polite");
    overlay.innerHTML = `
      <div class="working-dialog">
        <div class="working-spinner" aria-hidden="true"></div>
        <div class="working-copy">
          <h2 class="working-title">Working</h2>
          <p class="working-message">Starting local job.</p>
          <p class="working-elapsed">Elapsed time: 0:00</p>
          <div class="working-progress" aria-hidden="true">
            <span></span>
          </div>
        </div>
      </div>
    `;
    document.body.appendChild(overlay);
    return overlay;
  }

  function showWaitingScreen(form) {
    const currentOverlay = ensureOverlay();
    const title = form.dataset.waitTitle || "Working";
    const messages = messagesForForm(form);
    const titleNode = currentOverlay.querySelector(".working-title");
    const messageNode = currentOverlay.querySelector(".working-message");
    const elapsedNode = currentOverlay.querySelector(".working-elapsed");
    let messageIndex = 0;
    const startedAt = Date.now();

    clearInterval(elapsedTimer);
    clearInterval(messageTimer);

    titleNode.textContent = title;
    messageNode.textContent = messages[messageIndex];
    elapsedNode.textContent = "Elapsed time: 0:00";
    currentOverlay.classList.add("is-visible");
    document.body.classList.add("is-waiting");

    elapsedTimer = setInterval(() => {
      const elapsedSeconds = Math.floor((Date.now() - startedAt) / 1000);
      elapsedNode.textContent = `Elapsed time: ${formatElapsed(elapsedSeconds)}`;
    }, 1000);

    messageTimer = setInterval(() => {
      messageIndex = (messageIndex + 1) % messages.length;
      messageNode.textContent = messages[messageIndex];
    }, 9000);
  }

  function attachWaitingScreens() {
    const forms = document.querySelectorAll("form[data-wait-title]");
    for (const form of forms) {
      form.addEventListener("submit", (event) => {
        if (form.dataset.waitSubmitted === "true") {
          return;
        }

        event.preventDefault();
        form.dataset.waitSubmitted = "true";
        showWaitingScreen(form);

        for (const button of form.querySelectorAll("button[type='submit'], input[type='submit']")) {
          button.disabled = true;
        }

        window.setTimeout(() => {
          HTMLFormElement.prototype.submit.call(form);
        }, 80);
      });
    }
  }

  window.addEventListener("pageshow", () => {
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
