// Make banner divs clickable
document.querySelectorAll('.banner[data-href]').forEach(banner => {
  banner.addEventListener('click', () => {
    window.location.href = banner.getAttribute('data-href');
  });
});