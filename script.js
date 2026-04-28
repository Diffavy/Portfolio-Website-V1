const dialog = document.querySelectorAll('dialog');
const closeButton = document.querySelectorAll('.close-dialog');


closeButton.forEach(button => {
    button.addEventListener('click', () => {
        button.closest('dialog').close();
    });
});
 
dialog.forEach(d => {
    d.addEventListener("click", (e) => {
        if (e.target === d) {
            d.close()
        }
    })
})
// Add event listener to the project card to open the dialog
const projectCard = document.querySelectorAll('.project-card');

projectCard.forEach(card => {
    card.addEventListener('click', () => {
        card.nextElementSibling.showModal();
    });
});

// Theme switching
const themeSelector = document.getElementById('theme-selector');

if (localStorage.getItem('theme')) {
    const savedTheme = localStorage.getItem('theme')
    if (savedTheme === 'blueprint') {
        document.documentElement.removeAttribute('data-theme');
    } else {
        document.documentElement.setAttribute('data-theme', savedTheme)
    }
    themeSelector.value = savedTheme
}

 themeSelector.addEventListener('change', (event) => {
    const selectedTheme = event.target.value;
    localStorage.setItem('theme', selectedTheme)
    if (selectedTheme === 'blueprint') {
        document.documentElement.removeAttribute('data-theme');
    } else {
        document.documentElement.setAttribute('data-theme', selectedTheme);
    }
 })