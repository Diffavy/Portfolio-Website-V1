const dialog = document.getElementById('edutainment-dialog');
const closeButton = document.getElementById('close-dialog');

closeButton.addEventListener('click', () => {
    dialog.close();
});

// Add event listener to the project card to open the dialog
const projectCard = document.querySelector('.project-card');

projectCard.addEventListener('click', () => {
    dialog.showModal();
});