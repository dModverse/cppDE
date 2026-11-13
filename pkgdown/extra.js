// A code block of more than 15 lines folds into a details element; its summary says whether the
// block holds code, output (lines starting with #>) or both, and how many lines.
document.addEventListener("DOMContentLoaded", () => {
  const limit = 15;
  document.querySelectorAll("main pre").forEach((pre) => {
    const block = pre.closest("div.sourceCode") || pre;
    if (block.closest("details")) return;
    const lines = pre.textContent.replace(/\n+$/, "").split("\n");
    if (lines.length <= limit) return;
    const output = lines.filter((line) => line.startsWith("#>")).length;
    const kind = output === 0 ? "Code" : output === lines.length ? "Output" : "Code and output";
    const details = document.createElement("details");
    details.className = "fold";
    const summary = document.createElement("summary");
    summary.textContent = `${kind}, ${lines.length} lines`;
    details.append(summary);
    block.replaceWith(details);
    details.append(block);
  });
});
