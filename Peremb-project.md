src_dir: ./
output_dir: ./doc
docmark: >
predocmark: !
predocmark_alt: #
project: PEREMB
summary: A Subsystem DFT package able to combine arbitrary - periodic and non-periodic -  molecular systems
author: Christian Schwermann
author_description: PhD student in the [Doltsinis group](https://www.uni-muenster.de/Physik.FT/Forschung/agdoltsinis/index.html) at [WWU Muenster](https://www.uni-muenster.de)
email: c.schwermann@wwu.de
print_creation_date: true
year: 2019 Christian Schwermann <br/> ALL RIGHTS RESERVED
graph: true

# Usage #
PEREMB uses *densities* (as well as other molecular information, e.g. *density gradients*, *nuclear potential* or *nuclear positions*) 
of multiple subsystems as input and returns an effective *embedding potential* generated from inactive subsystems acting on an 
active subsystem as output.  

@note
This is a note
@endnote

> **Test**  
> blockquotes etc.

* lists
* nice
* code:
    call peremb( density, gradient, potential )
* doesnt seem to work
* equations!
    \begin{equation}
      \hat{H}\psi(\textbf{r},t) = i\hbar\frac{\partial}{\partial t}\psi(\textbf{r},t)  
    \end{equation}

Research group website: [Doltsinis group](https://www.uni-muenster.de/Physik.FT/Forschung/agdoltsinis/index.html)
