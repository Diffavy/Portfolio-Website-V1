const dialog = document.querySelectorAll('dialog');
const closeButton = document.querySelectorAll('.close-dialog');

// Close button functionality for dialogs
closeButton.forEach(button => {
    button.addEventListener('click', () => {
        button.closest('dialog').close();
        document.body.style.overflow = ''; // Re-enable background scrolling when dialog is closed
    });
});
 
// Close dialog when clicking outside of the content area
dialog.forEach(d => {
    d.addEventListener("click", (e) => {
        if (e.target === d) {
            d.close()
            document.body.style.overflow = ''; // Re-enable background scrolling when dialog is closed
        }
    })
})
// Add event listener to the project card to open the dialog
const projectCard = document.querySelectorAll('.project-card');

// Open the dialog when the project card is clicked
projectCard.forEach(card => {
    card.addEventListener('click', () => {
        const modal = card.nextElementSibling;
        modal.showModal();
        modal.scrollTop = 0; // Reset scroll position to top when opened
        modal.scrollIntoView({ behavior: 'smooth', block: 'center' }); // Scroll the dialog into view when opened
        setTimeout(() => {
            document.body.style.overflow = 'hidden'; // Disable background scrolling when dialog is open, delay to allow scrollIntoView to complete
        }, 500);
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


 // Cursor animation for hero section 

let i = 0;
const beforeCursor = "Harrun Muhiden ";
const afterCursor = " Software Engineer";
const divider = "<span class='cursor hero-divider'>|</span>";
const nameEl = document.querySelector('.hero-highlight');

const cursorEl = document.getElementById('type-title');

function typeTitle() {
    if (i < beforeCursor.length) {
        nameEl.innerHTML += beforeCursor[i]
        i++
        setTimeout(typeTitle, 100);
    } else if (i === beforeCursor.length) {
        nameEl.innerHTML += divider
        i++
        setTimeout(typeTitle, 200);
    } else if (i > (beforeCursor.length) && i <= (beforeCursor.length + afterCursor.length)) {
        nameEl.innerHTML += afterCursor[i - beforeCursor.length - 1]
        i++
        setTimeout(typeTitle, 100);
    } else {
        cursorEl.style.display = 'none'
        document.querySelector('.hero-divider').style.opacity = '1'
    }

}

typeTitle()

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

// Scroll animations
// For skill cards
const skillCards = document.querySelectorAll('.skill-card');

const skillCardObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('fade-in-scale');
        }})
    });

skillCards.forEach(card => {
    skillCardObserver.observe(card);
});


// For education cards
const educationCards = document.querySelectorAll('.education-card');

const educationCardObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('fade-in-down');
        }})
    });

educationCards.forEach(card => {
    educationCardObserver.observe(card);
});

// For contact section
const contactCards = document.querySelectorAll('.contact-card');

const contactCardObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('fade-in-scale');
        }})
    });

contactCards.forEach(card => {
    contactCardObserver.observe(card);
});

// For subheadings

const subheadings = document.querySelectorAll('.subheading');

const subheadingObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('subheading-entrance');
        }})
    });

subheadings.forEach(subheading => {
    subheadingObserver.observe(subheading);
});