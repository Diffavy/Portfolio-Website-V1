const dialog = document.querySelectorAll('dialog');
const closeButton = document.querySelectorAll('.close-dialog');

closeButton.forEach(button => {
    button.addEventListener('click', () => {
        button.closest('dialog').close();
    });
});
 
// Add event listener to the project card to open the dialog
const projectCard = document.querySelectorAll('.project-card');

projectCard.forEach(card => {
    card.addEventListener('click', () => {
        card.nextElementSibling.showModal();
    });
});