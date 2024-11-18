
Run scipti> main.m

Overview
We present a novel method for the automatic segmentation of the thoracic cavity and the detection of human lungs and the major thoracic organs, as a necessary pre-processing step for a subsequent deformable registration scheme. Our method is divided into two parts. In the first part, a coarse separation of the body into two sub-regions, one with the skeleton, fat and skin and another with the lungs, the heart and the rest soft tissues is achieved by approximating the rib cage with a c-spline curve. The second part employs curve-fitting methods in order to detect the boundaries of the lungs and the heart with better accuracy, based on the estimation of the first part.

Motivation

For a fast and efficient radiotherapy treatment planning of lung cancer, a successful registration procedure, which aligns the CT data and tracks accurately the thoracic volume during breathing cycle, is extremely essential for the radiologists.
Breathing motion changes dramatically the positions of the thoracic biological structures. During radiotherapy the position of the soft tissues and the tumor have to be defined as the x-rays beam direction has to be accurate. Registration procedure intends to capture the tissues deformations and displacements during breath cycle and to give information about the magnitude and the direction of these changes which will be used and evaluated from the clinician during planning process.
Currently deformable registration formulations fail to give accurate results when it is applied in the whole thoracic cavity. However, segmentation of the thoracic region into two subregions and then application of the registration algorithm in each subregion separately can give a solution.

In this project our aim was to segment the thoracic cavity as a pre-step for the registration process.
This work was a contribution to the image registration project for the tracking of lung tumors which is held at CMIC group -UCL (http://cmic.cs.ucl.ac.uk/).

Description
In the  first part we have a coarse distinction of the body into two sub-regions.
The one sub-region consists the outer thoracic part with the bones and the fat tissues and the other sub-region
(the inner thoracic part) with the lungs, the heart and some other tissues and structures like the aorta.
In detail, the "rough" separation into the two sub-regions is performed by estimating a curve which approximates the rib cage. A simple threshold selection technique is applied to segment the rib cage and c-spline interpolation is used to define the curve.


In the second part, we apply a coupled level set procedure in order to compensate for the "rough" separation of the first part and accurately extract  the whole thoracic cavity.
The level set framework employs a two phase scheme for the capture of the lungs and soft tissues (i.e. the heart and some other small tissues). The two-phase level set method intends to segment two objects (the lungs, the heart) and thus two zero level set functions should be defined, one for the lung and one for the heart. Initialization of the two zero level sets is achieved using the rib cage approximating curve.

For the initialization of the zero level function which will extract the heart we also have to define an initial small region inside the rib cage where the heart is positioned. For this purpose, a naive Bayes classifier is used to define which of the gray value pixels inside the estimated rib cage are most probable to belong to the heart. The pixels with the highest probability are used for the heart zero level function initialization. In the case of the lungs, the zero level function is initialized using the rib cage curve and the heart zero level set function.

Full description is here> https://www.researchgate.net/publication/304581790_Automatic_segmentation_of_the_thoracic_organs_for_image_registration_and_RT_planning
