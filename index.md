--- 
title: "Multi-omic data science with R/Bioconductor"
author: "University of Oulu & University of Turku"
date: "2022-04-23"
site: bookdown::bookdown_site
documentclass: book
bibliography: [packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: microbiome/course_2022_oulu
description: "Course material"
output:
  bookdown::gitbook
  bookdown::pdf_document2
always_allow_html: true  
classoption: oneside
geometry:
  - top=30mm
  - left=15mm
---


# Overview 

**Welcome to [Oulu Summer School, June 2022]()**

**Venue** University of Oulu. June 20-23, 2022.
  Organized together with [University of Turku](http://datascience.utu.fi), Finland.


<img src="https://user-images.githubusercontent.com/60338854/121848694-1072a480-ccf3-11eb-9af2-7fdefd8d1794.png" alt="ML4microbiome" width="50%"/>

<p style="font-size:12px">Figure source: Moreno-Indias _et al_. (2021) [Statistical and Machine Learning Techniques in Human Microbiome Studies: Contemporary Challenges and Solutions](https://doi.org/10.3389/fmicb.2021.635781). Frontiers in Microbiology 12:11.</p>


## Introduction

This course is based on data science with R, a popular open source
environment for scientific data analysis. It provides general
capabilities for the analysis, integration, and visualization of
multi-omic data from biomedical studies.

The course is based on [_miaverse_](https://microbiome.github.io) (mia = **MI**crobiome **A**nalysis), an
R/Bioconductor framework for microbiome data science. 

The data science framework consists of the following elements:

- efficient multi-omic data container 
- a package ecosystem, providing algorithmic data analysis methods
- demonstration data sets
- open documentation

The framework is explained in a greater depth in the online book
[Orchestrating Microbiome
Analysis](https://microbiome.github.io/OMA). The book is currently a
development version.


## Learning goals

The course aims to provide the basic understanding and skills for
biomedical data analysis with R and Bioconductor. The course provides
an overview of reproducible data analysis and reporting workflow in
multi-omic studies, with recent example data sets published microbiome
studies.

The participants become familiar with standard concepts and methods in
multi-omic data analysis, in particular in human microbiome research.
This includes better understanding of the specific statistical
challenges, practical hands-on experience with the commonly used
methods, and reproducible research with R.

After the course you will know how to approach new tasks in microbiome
data science by utilizing available documentation and R tools.

**Target audience** Advanced MSc and PhD students, Postdocs, and
  biomedical researchers who wish to develop their skills in
  scientific programming and biomedical data analysis.



## Schedule

The course takes place daily between 9am – 5pm, including coffee and
lunch breaks.

The course will be organized in a live format but the material will be
openly available online during and after the course.

A priority will given for local students from Oulu. Participants from
other higher education institutions are welcome to apply.

The mornings will start with lectures, and afternoons are mainly
dedicated to hands-on sessions that consist of practical tasks and
example data from recent research literature. Students will solve the
exercises based on available online examples and resources that are
pointed out in the study material. There is often more than one way to
solve a given task, and the teachers will be available for assistance.

Monday:

- Orientation
- Best practices in reproducible reporting and open science
- Hands-on: Introduction to R

Tuesday: 

- Key concepts and challenges in biomedical data analysis
- Hands-on: Biomedical data exploration 


Wednesday:

- Key concepts in biomedical data visualization
- Hands-on: Biomedical data visualization


Thursday:

- Advanced topics: common machine learning techniques
- Hands-on: multi-omic data integration and reproducible workflows


Friday:

- Student presentations
- Summary & Conclusions



## Material

The teaching material follows open online documentation created by the
course teachers, extending the online book Orchestrating Microbiome
Analysis (https://microbiome.github.io/OMA). We will teach generic
data analytical skills that are applicable to common data analysis
tasks encountered in modern omics research.

The training material walks through example workflows that go through
standard steps of biomedical data analysis covering data access,
exploration, analysis, visualization, reproducible reporting, and best
practices in open science. The teaching format allows adaptations
according to the student's learning speed.

**You can run the workflow by simply copy-pasting the
examples.** For further, advanced material, you can test and modify
further examples from the [OMA
book](https://microbiome.github.io/OMA), or try to apply the
techniques to your own data.

We expect that the students will install the necessary software in
advance. Online support will be available.



# Organizers

Jointly organized by:

- Health and Biosciences Doctoral Programme, University of Oulu Graduate School (HBS-DP)
- Department of Computing, University of Turku

Supported by:

- IT Center for Science (CSC), Finland

 
**Teachers**

  - [Leo Lahti](https://datascience.utu.fi) is the main teacher and
    Associate Professor in Data Science at the University of Turku,
    with specialization on biomedical data analysis.

  - Tuomas Borman is research assistant and one of the main developers
    of the open training material covered by the course.

  - Jenni Hekkala is a local PhD researchers from Oulu who will
    coordinate the local arrangements and contribute to course
    teaching.



## Acknowledgments

**Citation** "Introduction to microbiome data science (2021). URL: https://microbiome.github.io".

@oulu2022course


We thank all [miaverse developers and contributors](https://microbiome.github.io) who have contributed open resources that supported the development of the training material.

**Contact** [Leo Lahti](http://datascience.utu.fi), University of Turku, Finland

**License** All material is released under the open [CC BY-NC-SA 3.0 License](LICENSE) and available online during and after the course, following the
[recommendations on open teaching materials](https://avointiede.fi/fi/linjaukset-ja-aineistot/kotimaiset-linjaukset/oppimisen-ja-oppimateriaalien-avoimuuden-linjaus) of the national open science coordination in Finland**.

**Source code**

The source code of this repository is fully reproducible and contains
the Rmd files with executable code. All files can be rendered at one
go by running the file [main.R](main.R). You can check the file for
details on how to clone the repository and convert it into a gitbook,
although this is not necessary for the training.

- Source code (github): [miaverse teaching material](https://github.com/microbiome/course_2022_oulu)
- Course page (html): [miaverse teaching material](https://microbiome.github.io/course_2022_oulu/)

