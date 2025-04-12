# Genomics-WES-Project
## Day One: Introduction to Genetics, Genomics, and Precision Medicine

Welcome to Week One of our  series! Today’s focus is on establishing a solid foundation in genetics, genomics, and precision medicine—and understanding how sequencing data is generated and analyzed.

## Resources for Day One

Below are the key resources you need to review today:

### 1. Genetics, Genomics, and Precision Medicine
- **Resource:** [Genomics and Precision Medicine (YouTube)](https://www.youtube.com/watch?v=edXEzu9cGgg&t=1845s&pp=ygUfZ2Vub21pY3MgYW5kIHByZWNpc2lvbiBtZWRpY2luZQ%3D%3D)
- **Description:** This video introduces the basic concepts of genetics, genomics, and how precision medicine leverages genomic data to tailor medical treatment.

### 2. Illumina Introduction to Sequencing Data Analysis
- **Resource:** [Illumina Introduction to Sequencing Data Analysis (YouTube)](https://www.youtube.com/watch?v=XdoAnsDPMfA&pp=ygUgd2hvbGUgZXhvbWUgc2VxdWVuY2luZyBleHBsYWluZWQ%3D)
- **Description:** Learn about Illumina’s sequencing technology and an overview of the sequencing data analysis process.

### 3. Bioinformatics Workflow of Whole Exome Sequencing
- **Resource:** [GDC DNA Seq Variant Calling Pipeline Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
- **Description:** Explore a detailed bioinformatics workflow—from raw data through variant calling—for whole exome sequencing.

### 4. What is the Difference Between Genetics and Genomics?
- **Resource:** [Difference Between Genetics and Genomics (YouTube)](https://youtu.be/9p2q867t7Ow?si=mTme8SboDU0_w4p5)
- **Description:** This video clarifies the differences between genetics and genomics, which is crucial for understanding personalized medicine.

## Day 2: Today is a bit technical

**Task 0: Learn how to use Git and GitHub**
- **Instructions:**
1. Follow this tutorial on YT (https://youtu.be/tRZGeaHPoaw?si=vDTa13s5iRsiFpp6),
(https://youtu.be/b9fL9y_71d4?si=xwuyrjBrnE-I9rvK)
2. Post one question or insight in GitHub Discussions from what you learned on the first day on this repo** (https://github.com/ttobio/Genomics-WES-Project.git)

### Additional Resources - u01  
Here are some additional resources that helped me understand the concepts better:  
- [Markdown Guide](https://www.markdownguide.org/) - A comprehensive guide on Markdown syntax.  
- [GitHub Docs](https://docs.github.com/en) - Official GitHub documentation.  
- [Video Tutorial on GitHub & Markdown](https://youtu.be/fDkR0TDR9dI?si=vILUoDYfh2m4U4Cm) - A helpful video explaining GitHub and Markdown concepts.

**Task 1: Install Docker on WSL**
- **Objective:** Set up Docker on your Windows Subsystem for Linux.
- **Instructions:**  
  1. Follow this tutorial on YT (https://youtu.be/5RQbdMn04Oc?si=S1xmGl2okvKhWyJm) to install and configure Docker.
  2. Verify your installation by running `docker --version` in your terminal.

**Task 2: Install Bioinformatics Tools**
- **Tools:**  
  - **FastQC:** For quality control of sequencing data.
  - **BWA:** For aligning sequencing reads to a reference genome.
- **Instructions:**  
  1. Fastqc should be installed on your system but if not using conda `conda install bioconda::fastqc`
  2. Install BWA `conda install bioconda::bwa`
  2. Ensure that both tools run properly from the command line (e.g., try `fastqc --help` and `bwa`).

## Discussion and Feedback

- **Discussion:** Post any questions or installation issues in our GitHub Issues page.
- **Reflection:** Write a short reflection (2-3 sentences) about what you learned today and any challenges you encountered. Share this in the discussion forum.

## Day 3: Get a rest

**Task 1: Your first pull request**
- **Instructions**
1. Clone the repository on you local machine as a fork 
2. Create a new file in the assignments folder with your user id --> example u00.txt
3. Add and commit your changes and pull a request
- *You can type whatever you want in this file, just don't add confdential info*
- *This should be easy and it is covered in the YT videos provided*

**Task 2: Get to know markdown files**
- **Instructions**
1. [Get familiar with the markdown files](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
2. Apply this knowledge to the readme file (This file): maybe add resource that helped you and isn't here

- **Reflection:** Write a short reflection (2-3 sentences) about what you learned today and any challenges you encountered. Share this in the discussion forum.

## Additional Resources  - u04
- [Markdown Cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)  
- [Set up Git on macOS](https://www.atlassian.com/git/tutorials/install-git#mac-os-x)  
- [Basic Terminal Commands for macOS](https://www.macworld.com/article/671468/how-to-use-the-terminal-on-mac.html)

## Day 4: Know more about Docker and SHELL

**Task 1: Get used to Docker**
- **Instructions**
1. [Docker 101](https://youtu.be/eGz9DS-aIeY?si=KYJRg5YzS0JVQKpW)
2. [More of Docker](https://youtu.be/RqTEHSBrYFw?si=CSK5VMEb79hTFRzq)
3. Prepare a 3 minutes talk about Docker and how it will us in this project for next friday

**Task 2: Install GATK using Docker**
- **Instructions**
1. [Use this link to download GATK using Docker](https://hub.docker.com/r/broadinstitute/gatk/)
2. [More of GATK documentation](https://github.com/broadinstitute/gatk)
3. Verify the download and tell us about any tricky problems you faced on the discussion panel

**Task 3: Know more about shell**
- **Instructions**
1. [The Shell playlist](https://youtube.com/playlist?list=PLA86D04D6E0BFD2E0&si=bLIDog2l41W0w0Zv)


- **Reflection:** Write a short reflection (2-3 sentences) about what you learned today and any challenges you encountered. Share this in the discussion forum.

## Day 5: File formats
- **Instructions**
#### - Here are some YT videos that will help you understand the file formats we are going to work with 
1. [FASTA & FASTQ](https://youtu.be/xYEre9DtIqA?si=MsgBaiiCILqNry5-)
2. [VCF & gVCF](https://youtu.be/cBDZSGf9HRE?si=Oi5WkQCZrxhDuw40)
3. [SAM/BAM](https://youtu.be/vh6l9aPFUgU?si=F_T-ltv6TPh00Lc9)
4. [GFF3 & GTF](https://youtu.be/vh6l9aPFUgU?si=F_T-ltv6TPh00Lc9)


- **Reflection:** Post some interesting facts about the file formats that you learned today on the discussion panel


## Project-1: WES workflow 
- **Instructions**
### We will try to implement this workflow in GALAXY on our local machines [link](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html#data-preparation)
0. In the assignment folder, make a folder with the usernames of the group, for example --> u010522 and have this folder architecture:

├── workflow.md

├── output

├── plots

└── scripts --> might be scripts that helped you in the work flow (not single commands)
1. Download thr raw reads from this drive [link](https://drive.google.com/drive/folders/1iAThBHrCuCx-j6IPG60t74afDhmhzr41?usp=drive_link) 
#### Note: I would've loved to store the data on this repo but the limit for upload is 100Mb
#### There are 6 files for three poeple there --> mother, father, and proband
2. Check the tools in the tutorial and keep in mind the following:
  * Check if the tool is part of GATK and if not, does it have a replacement in GATK.
  * If not, install it your self on your local machine and include a link for the installation and the steps you followed.
3. Run each tool and explain the rationale behind each step and why you used a certain parameters for each tool. This should be included in your workflow.md
4. all the output files are stored in the output folder and sane idea for the plots. if you used any scripts from documentations  or a YT video, include the link the md file and put the script in the script folder.
5. The steps in the md file you will provide must be reproducable, meaning if someone has the tools installed and ran the same commands, it should produce the same output and the md file should contain all the details for someone to be able to both follow it and understand the rationale behind it.


- **Reminder:** If you had any problems, search the internet and if you didn't find a solution, feel free to post it in the issues or on the WhatsApp group.
