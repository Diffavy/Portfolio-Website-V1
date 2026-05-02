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
        const modal = card.nextElementSibling;
        modal.showModal();
        modal.scrollTop = 0;
    });
});

// Theme switching
const themeSelector = document.getElementById('theme-selector');

if (localStorage.getItem('theme')) {
    const savedTheme = localStorage.getItem('theme')
    if (savedTheme) {
        themeSelector.value = savedTheme
    }
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

 // Set the current year in the footer
 const yearElement = document.getElementById('current-year');
 yearElement.textContent = new Date().getFullYear();

// Set the age in the about me section
const dob = new Date('2003-04-12')
const today = new Date()
let age = today.getFullYear() - dob.getFullYear()

const hasHadBirthdayThisYear = 
    today.getMonth() > dob.getMonth() || 
    (today.getMonth() === dob.getMonth() && today.getDate() >= dob.getDate())

if (!hasHadBirthdayThisYear) {
    age--
}

document.getElementById('age').textContent = age