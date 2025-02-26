# Latest Papers on Astrophysics (2025-02-26)

## CLBlast: A Tuned OpenCL BLAS Library
**Published:** 2017-05-12
**Authors:** Cedric Nugteren
**Abstract:** This work introduces CLBlast, an open-source BLAS library providing optimized
OpenCL routines to accelerate dense linear algebra for a wide variety of
devices. It is targeted at machine learning and HPC applications and thus
provides a fast matrix-multiplication routine (GEMM) to accelerate the core of
many applications (e.g. deep learning, iterative solvers, astrophysics,
computational fluid dynamics, quantum chemistry). CLBlast has five main
advantages over other OpenCL BLAS libraries: 1) it is optimized for and tested
on a large variety of OpenCL devices including less commonly used devices such
as embedded and low-power GPUs, 2) it can be explicitly tuned for specific
problem-sizes on specific hardware platforms, 3) it can perform operations in
half-precision floating-point FP16 saving bandwidth, time and energy, 4) it has
an optional CUDA back-end, 5) and it can combine multiple operations in a
single batched routine, accelerating smaller problems significantly. This paper
describes the library and demonstrates the advantages of CLBlast experimentally
for different use-cases on a wide variety of OpenCL hardware.
**URL:** No URL available

## Advice from the Oracle: Really Intelligent Information Retrieval
**Published:** 2018-01-02
**Authors:** Michael J. Kurtz
**Abstract:** What is "intelligent" information retrieval? Essentially this is asking what
is intelligence, in this article I will attempt to show some of the aspects of
human intelligence, as related to information retrieval. I will do this by the
device of a semi-imaginary Oracle. Every Observatory has an oracle, someone who
is a distinguished scientist, has great administrative responsibilities, acts
as mentor to a number of less senior people, and as trusted advisor to even the
most accomplished scientists, and knows essentially everyone in the field. In
an appendix I will present a brief summary of the Statistical Factor Space
method for text indexing and retrieval, and indicate how it will be used in the
Astrophysics Data System Abstract Service. 2018 Keywords: Personal Digital
Assistant; Supervised Topic Models
**URL:** No URL available

## Deep Learning for Real-time Gravitational Wave Detection and Parameter Estimation with LIGO Data
**Published:** 2017-11-21
**Authors:** E. A. Huerta, Daniel George
**Abstract:** The recent Nobel-prize-winning detections of gravitational waves from merging
black holes and the subsequent detection of the collision of two neutron stars
in coincidence with electromagnetic observations have inaugurated a new era of
multimessenger astrophysics. To enhance the scope of this emergent science, we
proposed the use of deep convolutional neural networks for the detection and
characterization of gravitational wave signals in real-time. This method, Deep
Filtering, was initially demonstrated using simulated LIGO noise. In this
article, we present the extension of Deep Filtering using real data from the
first observing run of LIGO, for both detection and parameter estimation of
gravitational waves from binary black hole mergers with continuous data streams
from multiple LIGO detectors. We show for the first time that machine learning
can detect and estimate the true parameters of a real GW event observed by
LIGO. Our comparisons show that Deep Filtering is far more computationally
efficient than matched-filtering, while retaining similar sensitivity and lower
errors, allowing real-time processing of weak time-series signals in
non-stationary non-Gaussian noise, with minimal resources, and also enables the
detection of new classes of gravitational wave sources that may go unnoticed
with existing detection algorithms. This approach is uniquely suited to enable
coincident detection campaigns of gravitational waves and their multimessenger
counterparts in real-time.
**URL:** No URL available

## Denoising Gravitational Waves using Deep Learning with Recurrent Denoising Autoencoders
**Published:** 2017-11-27
**Authors:** Hongyu Shen, Zhizhen Zhao, E. A. Huerta, Daniel George
**Abstract:** Gravitational wave astronomy is a rapidly growing field of modern
astrophysics, with observations being made frequently by the LIGO detectors.
Gravitational wave signals are often extremely weak and the data from the
detectors, such as LIGO, is contaminated with non-Gaussian and non-stationary
noise, often containing transient disturbances which can obscure real signals.
Traditional denoising methods, such as principal component analysis and
dictionary learning, are not optimal for dealing with this non-Gaussian noise,
especially for low signal-to-noise ratio gravitational wave signals.
Furthermore, these methods are computationally expensive on large datasets. To
overcome these issues, we apply state-of-the-art signal processing techniques,
based on recent groundbreaking advancements in deep learning, to denoise
gravitational wave signals embedded either in Gaussian noise or in real LIGO
noise. We introduce SMTDAE, a Staired Multi-Timestep Denoising Autoencoder,
based on sequence-to-sequence bi-directional Long-Short-Term-Memory recurrent
neural networks. We demonstrate the advantages of using our unsupervised deep
learning approach and show that, after training only using simulated Gaussian
noise, SMTDAE achieves superior recovery performance for gravitational wave
signals embedded in real non-Gaussian LIGO noise.
**URL:** No URL available

## Deep Neural Networks to Enable Real-time Multimessenger Astrophysics
**Published:** 2016-12-30
**Authors:** E. A. Huerta, Daniel George
**Abstract:** Gravitational wave astronomy has set in motion a scientific revolution. To
further enhance the science reach of this emergent field, there is a pressing
need to increase the depth and speed of the gravitational wave algorithms that
have enabled these groundbreaking discoveries. To contribute to this effort, we
introduce Deep Filtering, a new highly scalable method for end-to-end
time-series signal processing, based on a system of two deep convolutional
neural networks, which we designed for classification and regression to rapidly
detect and estimate parameters of signals in highly noisy time-series data
streams. We demonstrate a novel training scheme with gradually increasing noise
levels, and a transfer learning procedure between the two networks. We showcase
the application of this method for the detection and parameter estimation of
gravitational waves from binary black hole mergers. Our results indicate that
Deep Filtering significantly outperforms conventional machine learning
techniques, achieves similar performance compared to matched-filtering while
being several orders of magnitude faster thus allowing real-time processing of
raw big data with minimal resources. More importantly, Deep Filtering extends
the range of gravitational wave signals that can be detected with ground-based
gravitational wave detectors. This framework leverages recent advances in
artificial intelligence algorithms and emerging hardware architectures, such as
deep-learning-optimized GPUs, to facilitate real-time searches of gravitational
wave sources and their electromagnetic and astro-particle counterparts.
**URL:** No URL available

## Deep Learning for Real-time Gravitational Wave Detection and Parameter Estimation: Results with Advanced LIGO Data
**Published:** 2017-11-08
**Authors:** E. A. Huerta, Daniel George
**Abstract:** The recent Nobel-prize-winning detections of gravitational waves from merging
black holes and the subsequent detection of the collision of two neutron stars
in coincidence with electromagnetic observations have inaugurated a new era of
multimessenger astrophysics. To enhance the scope of this emergent field of
science, we pioneered the use of deep learning with convolutional neural
networks, that take time-series inputs, for rapid detection and
characterization of gravitational wave signals. This approach, Deep Filtering,
was initially demonstrated using simulated LIGO noise. In this article, we
present the extension of Deep Filtering using real data from LIGO, for both
detection and parameter estimation of gravitational waves from binary black
hole mergers using continuous data streams from multiple LIGO detectors. We
demonstrate for the first time that machine learning can detect and estimate
the true parameters of real events observed by LIGO. Our results show that Deep
Filtering achieves similar sensitivities and lower errors compared to
matched-filtering while being far more computationally efficient and more
resilient to glitches, allowing real-time processing of weak time-series
signals in non-stationary non-Gaussian noise with minimal resources, and also
enables the detection of new classes of gravitational wave sources that may go
unnoticed with existing detection algorithms. This unified framework for data
analysis is ideally suited to enable coincident detection campaigns of
gravitational waves and their multimessenger counterparts in real-time.
**URL:** No URL available

## Convolutional Networks for Spherical Signals
**Published:** 2017-09-14
**Authors:** Jonas Köhler, Mario Geiger, Taco Cohen, Max Welling
**Abstract:** The success of convolutional networks in learning problems involving planar
signals such as images is due to their ability to exploit the translation
symmetry of the data distribution through weight sharing. Many areas of science
and egineering deal with signals with other symmetries, such as rotation
invariant data on the sphere. Examples include climate and weather science,
astrophysics, and chemistry. In this paper we present spherical convolutional
networks. These networks use convolutions on the sphere and rotation group,
which results in rotational weight sharing and rotation equivariance. Using a
synthetic spherical MNIST dataset, we show that spherical convolutional
networks are very effective at dealing with rotationally invariant
classification problems.
**URL:** No URL available

## Accelerating Approximate Bayesian Computation with Quantile Regression: Application to Cosmological Redshift Distributions
**Published:** 2017-07-24
**Authors:** Alexandre Réfrégier, Jörg Herbel, Tomasz Kacprzak, Adam Amara
**Abstract:** Approximate Bayesian Computation (ABC) is a method to obtain a posterior
distribution without a likelihood function, using simulations and a set of
distance metrics. For that reason, it has recently been gaining popularity as
an analysis tool in cosmology and astrophysics. Its drawback, however, is a
slow convergence rate. We propose a novel method, which we call qABC, to
accelerate ABC with Quantile Regression. In this method, we create a model of
quantiles of distance measure as a function of input parameters. This model is
trained on a small number of simulations and estimates which regions of the
prior space are likely to be accepted into the posterior. Other regions are
then immediately rejected. This procedure is then repeated as more simulations
are available. We apply it to the practical problem of estimation of redshift
distribution of cosmological samples, using forward modelling developed in
previous work. The qABC method converges to nearly same posterior as the basic
ABC. It uses, however, only 20\% of the number of simulations compared to basic
ABC, achieving a fivefold gain in execution time for our problem. For other
problems the acceleration rate may vary; it depends on how close the prior is
to the final posterior. We discuss possible improvements and extensions to this
method.
**URL:** No URL available

## By-passing the Kohn-Sham equations with machine learning
**Published:** 2016-09-09
**Authors:** Klaus-Robert Müller, Kieron Burke, Mark E. Tuckerman, Li Li, Leslie Vogt, Felix Brockherde
**Abstract:** Last year, at least 30,000 scientific papers used the Kohn-Sham scheme of
density functional theory to solve electronic structure problems in a wide
variety of scientific fields, ranging from materials science to biochemistry to
astrophysics. Machine learning holds the promise of learning the kinetic energy
functional via examples, by-passing the need to solve the Kohn-Sham equations.
This should yield substantial savings in computer time, allowing either larger
systems or longer time-scales to be tackled, but attempts to machine-learn this
functional have been limited by the need to find its derivative. The present
work overcomes this difficulty by directly learning the density-potential and
energy-density maps for test systems and various molecules. Both improved
accuracy and lower computational cost with this method are demonstrated by
reproducing DFT energies for a range of molecular geometries generated during
molecular dynamics simulations. Moreover, the methodology could be applied
directly to quantum chemical calculations, allowing construction of density
functionals of quantum-chemical accuracy.
**URL:** No URL available

## Big Universe, Big Data: Machine Learning and Image Analysis for Astronomy
**Published:** 2017-04-15
**Authors:** Kim Steenstrup Pedersen, Fabian Gieseke, Kristoffer Stensbo-Smidt, Jan Kremer, Christian Igel
**Abstract:** Astrophysics and cosmology are rich with data. The advent of wide-area
digital cameras on large aperture telescopes has led to ever more ambitious
surveys of the sky. Data volumes of entire surveys a decade ago can now be
acquired in a single night and real-time analysis is often desired. Thus,
modern astronomy requires big data know-how, in particular it demands highly
efficient machine learning and image analysis algorithms. But scalability is
not the only challenge: Astronomy applications touch several current machine
learning research questions, such as learning from biased data and dealing with
label and measurement noise. We argue that this makes astronomy a great domain
for computer science research, as it pushes the boundaries of data analysis. In
the following, we will present this exciting application area for data
scientists. We will focus on exemplary results, discuss main challenges, and
highlight some recent methodological advancements in machine learning and image
analysis triggered by astronomical applications.
**URL:** No URL available

## Phylogenetic Tools in Astrophysics
**Published:** 2017-03-01
**Authors:** Didier Fraix-Burnet
**Abstract:** Multivariate clustering in astrophysics is a recent development justified by
the bigger and bigger surveys of the sky. The phylogenetic approach is probably
the most unexpected technique that has appeared for the unsupervised
classification of galaxies, stellar populations or globular clusters. On one
side, this is a somewhat natural way of classifying astrophysical entities
which are all evolving objects. On the other side, several conceptual and
practical difficulties arize, such as the hierarchical representation of the
astrophysical diversity, the continuous nature of the parameters, and the
adequation of the result to the usual practice for the physical interpretation.
Most of these have now been solved through the studies of limited samples of
stellar clusters and galaxies. Up to now, only the Maximum Parsimony
(cladistics) has been used since it is the simplest and most general
phylogenetic technique. Probabilistic and network approaches are obvious
extensions that should be explored in the future.
**URL:** No URL available

## Selective De-noising of Sparse-Coloured Images
**Published:** 2016-10-29
**Authors:** Arjun Chaudhuri
**Abstract:** Since time immemorial, noise has been a constant source of disturbance to the
various entities known to mankind. Noise models of different kinds have been
developed to study noise in more detailed fashion over the years. Image
processing, particularly, has extensively implemented several algorithms to
reduce noise in photographs and pictorial documents to alleviate the effect of
noise. Images with sparse colours-lesser number of distinct colours in them-are
common nowadays, especially in astronomy and astrophysics where black and white
colours form the main components. Additive noise of Gaussian type is the most
common form of noise to be studied and analysed in majority of communication
channels, namely-satellite links, mobile base station to local cellular tower
communication channel,et. al. Most of the time, we encounter images from
astronomical sources being distorted with noise maximally as they travel long
distance from telescopes in outer space to Earth. Considering Additive White
Gaussian Noise(AWGN) to be the common noise in these long distance channels,
this paper provides an insight and an algorithmic approach to pixel-specific
de-noising of sparse-coloured images affected by AWGN. The paper concludes with
some essential future avenues and applications of this de-noising method in
industry and academia.
**URL:** No URL available

## Clustering with phylogenetic tools in astrophysics
**Published:** 2016-06-01
**Authors:** Didier Fraix-Burnet
**Abstract:** Phylogenetic approaches are finding more and more applications outside the
field of biology. Astrophysics is no exception since an overwhelming amount of
multivariate data has appeared in the last twenty years or so. In particular,
the diversification of galaxies throughout the evolution of the Universe quite
naturally invokes phylogenetic approaches. We have demonstrated that Maximum
Parsimony brings useful astrophysical results, and we now proceed toward the
analyses of large datasets for galaxies. In this talk I present how we solve
the major difficulties for this goal: the choice of the parameters, their
discretization, and the analysis of a high number of objects with an
unsupervised NP-hard classification technique like cladistics. 1. Introduction
How do the galaxy form, and when? How did the galaxy evolve and transform
themselves to create the diversity we observe? What are the progenitors to
present-day galaxies? To answer these big questions, observations throughout
the Universe and the physical modelisation are obvious tools. But between
these, there is a key process, without which it would be impossible to extract
some digestible information from the complexity of these systems. This is
classification. One century ago, galaxies were discovered by Hubble. From
images obtained in the visible range of wavelengths, he synthetised his
observations through the usual process: classification. With only one parameter
(the shape) that is qualitative and determined with the eye, he found four
categories: ellipticals, spirals, barred spirals and irregulars. This is the
famous Hubble classification. He later hypothetized relationships between these
classes, building the Hubble Tuning Fork. The Hubble classification has been
refined, notably by de Vaucouleurs, and is still used as the only global
classification of galaxies. Even though the physical relationships proposed by
Hubble are not retained any more, the Hubble Tuning Fork is nearly always used
to represent the classification of the galaxy diversity under its new name the
Hubble sequence (e.g. Delgado-Serrano, 2012). Its success is impressive and can
be understood by its simplicity, even its beauty, and by the many correlations
found between the morphology of galaxies and their other properties. And one
must admit that there is no alternative up to now, even though both the Hubble
classification and diagram have been recognised to be unsatisfactory. Among the
most obvious flaws of this classification, one must mention its monovariate,
qualitative, subjective and old-fashioned nature, as well as the difficulty to
characterise the morphology of distant galaxies. The first two most significant
multivariate studies were by Watanabe et al. (1985) and Whitmore (1984). Since
the year 2005, the number of studies attempting to go beyond the Hubble
classification has increased largely. Why, despite of this, the Hubble
classification and its sequence are still alive and no alternative have yet
emerged (Sandage, 2005)? My feeling is that the results of the multivariate
analyses are not easily integrated into a one-century old practice of modeling
the observations. In addition, extragalactic objects like galaxies, stellar
clusters or stars do evolve. Astronomy now provides data on very distant
objects, raising the question of the relationships between those and our
present day nearby galaxies. Clearly, this is a phylogenetic problem.
Astrocladistics 1 aims at exploring the use of phylogenetic tools in
astrophysics (Fraix-Burnet et al., 2006a,b). We have proved that Maximum
Parsimony (or cladistics) can be applied in astrophysics and provides a new
exploration tool of the data (Fraix-Burnet et al., 2009, 2012, Cardone \&
Fraix-Burnet, 2013). As far as the classification of galaxies is concerned, a
larger number of objects must now be analysed. In this paper, I
**URL:** No URL available

## Probabilistic Numerics and Uncertainty in Computations
**Published:** 2015-06-03
**Authors:** Michael A. Osborne, Philipp Hennig, Mark Girolami
**Abstract:** We deliver a call to arms for probabilistic numerical methods: algorithms for
numerical tasks, including linear algebra, integration, optimization and
solving differential equations, that return uncertainties in their
calculations. Such uncertainties, arising from the loss of precision induced by
numerical calculation with limited time or hardware, are important for much
contemporary science and industry. Within applications such as climate science
and astrophysics, the need to make decisions on the basis of computations with
large and complex data has led to a renewed focus on the management of
numerical uncertainty. We describe how several seminal classic numerical
methods can be interpreted naturally as probabilistic inference. We then show
that the probabilistic view suggests new algorithms that can flexibly be
adapted to suit application specifics, while delivering improved empirical
performance. We provide concrete illustrations of the benefits of probabilistic
numeric algorithms on real scientific problems from astrometry and astronomical
imaging, while highlighting open problems with these new algorithms. Finally,
we describe how probabilistic numerical methods provide a coherent framework
for identifying the uncertainty in calculations performed with a combination of
numerical algorithms (e.g. both numerical optimisers and differential equation
solvers), potentially allowing the diagnosis (and control) of error sources in
computations.
**URL:** No URL available

## Data Science and Ebola
**Published:** 2015-04-11
**Authors:** Aske Plaat
**Abstract:** Data Science---Today, everybody and everything produces data. People produce
large amounts of data in social networks and in commercial transactions.
Medical, corporate, and government databases continue to grow. Sensors continue
to get cheaper and are increasingly connected, creating an Internet of Things,
and generating even more data. In every discipline, large, diverse, and rich
data sets are emerging, from astrophysics, to the life sciences, to the
behavioral sciences, to finance and commerce, to the humanities and to the
arts. In every discipline people want to organize, analyze, optimize and
understand their data to answer questions and to deepen insights. The science
that is transforming this ocean of data into a sea of knowledge is called data
science. This lecture will discuss how data science has changed the way in
which one of the most visible challenges to public health is handled, the 2014
Ebola outbreak in West Africa.
**URL:** No URL available

## Machine Learning Etudes in Astrophysics: Selection Functions for Mock Cluster Catalogs
**Published:** 2014-09-04
**Authors:** J. Richard Bond, Amir Hajian, Marcelo Alvarez
**Abstract:** Making mock simulated catalogs is an important component of astrophysical
data analysis. Selection criteria for observed astronomical objects are often
too complicated to be derived from first principles. However the existence of
an observed group of objects is a well-suited problem for machine learning
classification. In this paper we use one-class classifiers to learn the
properties of an observed catalog of clusters of galaxies from ROSAT and to
pick clusters from mock simulations that resemble the observed ROSAT catalog.
We show how this method can be used to study the cross-correlations of thermal
Sunya'ev-Zeldovich signals with number density maps of X-ray selected cluster
catalogs. The method reduces the bias due to hand-tuning the selection function
and is readily scalable to large catalogs with a high-dimensional space of
astrophysical features.
**URL:** No URL available

## Spectral classification using convolutional neural networks
**Published:** 2014-12-29
**Authors:** Pavel Hála
**Abstract:** There is a great need for accurate and autonomous spectral classification
methods in astrophysics. This thesis is about training a convolutional neural
network (ConvNet) to recognize an object class (quasar, star or galaxy) from
one-dimension spectra only. Author developed several scripts and C programs for
datasets preparation, preprocessing and postprocessing of the data. EBLearn
library (developed by Pierre Sermanet and Yann LeCun) was used to create
ConvNets. Application on dataset of more than 60000 spectra yielded success
rate of nearly 95%. This thesis conclusively proved great potential of
convolutional neural networks and deep learning methods in astrophysics.
**URL:** No URL available

## Sparsity and adaptivity for the blind separation of partially correlated sources
**Published:** 2014-12-09
**Authors:** Jean-Luc Starck, Jeremy Rapin, Jerome Bobin, Anthony Larue
**Abstract:** Blind source separation (BSS) is a very popular technique to analyze
multichannel data. In this context, the data are modeled as the linear
combination of sources to be retrieved. For that purpose, standard BSS methods
all rely on some discrimination principle, whether it is statistical
independence or morphological diversity, to distinguish between the sources.
However, dealing with real-world data reveals that such assumptions are rarely
valid in practice: the signals of interest are more likely partially
correlated, which generally hampers the performances of standard BSS methods.
In this article, we introduce a novel sparsity-enforcing BSS method coined
Adaptive Morphological Component Analysis (AMCA), which is designed to retrieve
sparse and partially correlated sources. More precisely, it makes profit of an
adaptive re-weighting scheme to favor/penalize samples based on their level of
correlation. Extensive numerical experiments have been carried out which show
that the proposed method is robust to the partial correlation of sources while
standard BSS techniques fail. The AMCA algorithm is evaluated in the field of
astrophysics for the separation of physical components from microwave data.
**URL:** No URL available

## NMF with Sparse Regularizations in Transformed Domains
**Published:** 2014-07-29
**Authors:** Jean-Luc Starck, Jérôme Bobin, Jérémy Rapin, Anthony Larue
**Abstract:** Non-negative blind source separation (non-negative BSS), which is also
referred to as non-negative matrix factorization (NMF), is a very active field
in domains as different as astrophysics, audio processing or biomedical signal
processing. In this context, the efficient retrieval of the sources requires
the use of signal priors such as sparsity. If NMF has now been well studied
with sparse constraints in the direct domain, only very few algorithms can
encompass non-negativity together with sparsity in a transformed domain since
simultaneously dealing with two priors in two different domains is challenging.
In this article, we show how a sparse NMF algorithm coined non-negative
generalized morphological component analysis (nGMCA) can be extended to impose
non-negativity in the direct domain along with sparsity in a transformed
domain, with both analysis and synthesis formulations. To our knowledge, this
work presents the first comparison of analysis and synthesis priors ---as well
as their reweighted versions--- in the context of blind source separation.
Comparisons with state-of-the-art NMF algorithms on realistic data show the
efficiency as well as the robustness of the proposed algorithms.
**URL:** No URL available

## Pattern recognition issues on anisotropic smoothed particle hydrodynamics
**Published:** 2013-08-06
**Authors:** Eraldo Pereira Marinho
**Abstract:** This is a preliminary theoretical discussion on the computational
requirements of the state of the art smoothed particle hydrodynamics (SPH) from
the optics of pattern recognition and artificial intelligence. It is pointed
out in the present paper that, when including anisotropy detection to improve
resolution on shock layer, SPH is a very peculiar case of unsupervised machine
learning. On the other hand, the free particle nature of SPH opens an
opportunity for artificial intelligence to study particles as agents acting in
a collaborative framework in which the timed outcomes of a fluid simulation
forms a large knowledge base, which might be very attractive in computational
astrophysics phenomenological problems like self-propagating star formation.
**URL:** No URL available

## Infinite Shift-invariant Grouped Multi-task Learning for Gaussian Processes
**Published:** 2012-03-05
**Authors:** Roni Khardon, Yuyang Wang, Pavlos Protopapas
**Abstract:** Multi-task learning leverages shared information among data sets to improve
the learning performance of individual tasks. The paper applies this framework
for data where each task is a phase-shifted periodic time series. In
particular, we develop a novel Bayesian nonparametric model capturing a mixture
of Gaussian processes where each task is a sum of a group-specific function and
a component capturing individual variation, in addition to each task being
phase shifted. We develop an efficient \textsc{em} algorithm to learn the
parameters of the model. As a special case we obtain the Gaussian mixture model
and \textsc{em} algorithm for phased-shifted periodic time series. Furthermore,
we extend the proposed model by using a Dirichlet Process prior and thereby
leading to an infinite mixture model that is capable of doing automatic model
selection. A Variational Bayesian approach is developed for inference in this
model. Experiments in regression, classification and class discovery
demonstrate the performance of the proposed models using both synthetic data
and real-world time series data from astrophysics. Our methods are particularly
useful when the time series are sparsely and non-synchronously sampled.
**URL:** No URL available

## Fast Automated Analysis of Strong Gravitational Lenses with Convolutional Neural Networks
**Published:** 2017-08-29
**Authors:** Laurence Perreault Levasseur, Yashar D. Hezaveh, Philip J. Marshall
**Abstract:** Quantifying image distortions caused by strong gravitational lensing and
estimating the corresponding matter distribution in lensing galaxies has been
primarily performed by maximum likelihood modeling of observations. This is
typically a time and resource-consuming procedure, requiring sophisticated
lensing codes, several data preparation steps, and finding the maximum
likelihood model parameters in a computationally expensive process with
downhill optimizers. Accurate analysis of a single lens can take up to a few
weeks and requires the attention of dedicated experts. Tens of thousands of new
lenses are expected to be discovered with the upcoming generation of ground and
space surveys, the analysis of which can be a challenging task. Here we report
the use of deep convolutional neural networks to accurately estimate lensing
parameters in an extremely fast and automated way, circumventing the
difficulties faced by maximum likelihood methods. We also show that lens
removal can be made fast and automated using Independent Component Analysis of
multi-filter imaging data. Our networks can recover the parameters of the
Singular Isothermal Ellipsoid density profile, commonly used to model strong
lensing systems, with an accuracy comparable to the uncertainties of
sophisticated models, but about ten million times faster: 100 systems in
approximately 1s on a single graphics processing unit. These networks can
provide a way for non-experts to obtain lensing parameter estimates for large
samples of data. Our results suggest that neural networks can be a powerful and
fast alternative to maximum likelihood procedures commonly used in
astrophysics, radically transforming the traditional methods of data reduction
and analysis.
**URL:** No URL available

## Topological Data Analysis Made Easy with the Topology ToolKit
**Published:** 2018-06-21
**Authors:** Guillaume Favelier, Julien Tierny, Daisuke Sakurai, Charles Gueunet, Will Usher, Qi Wu, Joshua Levine, Maxime Soler, Attila Gyulassy, Julien Kitware, Jonas Lukasczyk
**Abstract:** This tutorial presents topological methods for the analysis and visualization
of scientific data from a user's perspective, with the Topology ToolKit (TTK),
a recently released open-source library for topological data analysis.
Topological methods have gained considerably in popularity and maturity over
the last twenty years and success stories of established methods have been
documented in a wide range of applications (combustion, chemistry,
astrophysics, material sciences, etc.) with both acquired and simulated data,
in both post-hoc and in-situ contexts. While reference textbooks have been
published on the topic, no tutorial at IEEE VIS has covered this area in recent
years, and never at a software level and from a user's point-of-view. This
tutorial fills this gap by providing a beginner's introduction to topological
methods for practitioners, researchers, students, and lecturers. In particular,
instead of focusing on theoretical aspects and algorithmic details, this
tutorial focuses on how topological methods can be useful in practice for
concrete data analysis tasks such as segmentation, feature extraction or
tracking. The tutorial describes in detail how to achieve these tasks with TTK.
First, after an introduction to topological methods and their application in
data analysis, a brief overview of TTK's main entry point for end users, namely
ParaView, will be presented. Second, an overview of TTK's main features will be
given. A running example will be described in detail, showcasing how to access
TTK's features via ParaView, Python, VTK/C++, and C++. Third, hands-on sessions
will concretely show how to use TTK in ParaView for multiple, representative
data analysis tasks. Fourth, the usage of TTK will be presented for developers,
in particular by describing several examples of visualization and data analysis
projects that were built on top of TTK. Finally, some feedback regarding the
usage of TTK as a teaching platform for topological analysis will be given.
Presenters of this tutorial include experts in topological methods, core
authors of TTK as well as active users, coming from academia, labs, or
industry. A large part of the tutorial will be dedicated to hands-on exercises
and a rich material package (including TTK pre-installs in virtual machines,
code, data, demos, video tutorials, etc.) will be provided to the participants.
This tutorial mostly targets students, practitioners and researchers who are
not experts in topological methods but who are interested in using them in
their daily tasks. We also target researchers already familiar to topological
methods and who are interested in using or contributing to TTK.
**URL:** No URL available

## swordfish: Efficient Forecasting of New Physics Searches without Monte Carlo
**Published:** 2017-12-14
**Authors:** Christoph Weniger, Thomas D. P. Edwards
**Abstract:** We introduce swordfish, a Monte-Carlo-free Python package to predict expected
exclusion limits, the discovery reach and expected confidence contours for a
large class of experiments relevant for particle- and astrophysics. The tool is
applicable to any counting experiment, supports general correlated background
uncertainties, and gives exact results in both the signal- and
systematics-limited regimes. Instead of time-intensive Monte Carlo simulations
and likelihood maximization, it internally utilizes new approximation methods
that are built on information geometry. Out of the box, swordfish provides
straightforward methods for accurately deriving many of the common sensitivity
measures. In addition, it allows one to examine experimental abilities in great
detail by employing the notion of information flux. This new concept
generalizes signal-to-noise ratios to situations where background uncertainties
and component mixing cannot be neglected. The user interface of swordfish is
designed with ease-of-use in mind, which we demonstrate by providing typical
examples from indirect and direct dark matter searches as jupyter notebooks.
**URL:** No URL available

## emcee: The MCMC Hammer
**Published:** 2012-02-16
**Authors:** Daniel Foreman-Mackey, Jonathan Goodman, Dustin Lang, David W. Hogg
**Abstract:** We introduce a stable, well tested Python implementation of the
affine-invariant ensemble sampler for Markov chain Monte Carlo (MCMC) proposed
by Goodman & Weare (2010). The code is open source and has already been used in
several published projects in the astrophysics literature. The algorithm behind
emcee has several advantages over traditional MCMC sampling methods and it has
excellent performance as measured by the autocorrelation time (or function
calls per independent sample). One major advantage of the algorithm is that it
requires hand-tuning of only 1 or 2 parameters compared to $\sim N^2$ for a
traditional algorithm in an N-dimensional parameter space. In this document, we
describe the algorithm and the details of our implementation and API.
Exploiting the parallelism of the ensemble method, emcee permits any user to
take advantage of multiple CPU cores without extra effort. The code is
available online at http://dan.iel.fm/emcee under the MIT License.
**URL:** No URL available

## Double Compact Objects III: Gravitational Wave Detection Rates
**Published:** 2014-05-27
**Authors:** K. Belczynski, R. O'Shaughnessy, I. Mandel, C. Fryer, E. Berti, T. Bulik, M. Dominik, F. Pannarale, D. Holz
**Abstract:** The unprecedented range of second-generation gravitational-wave (GW)
observatories calls for refining the predictions of potential sources and
detection rates. The coalescence of double compact objects (DCOs)---i.e.,
neutron star-neutron star (NS-NS), black hole-neutron star (BH-NS), and black
hole-black hole (BH-BH) binary systems---is the most promising source of GWs
for these detectors. We compute detection rates of coalescing DCOs in
second-generation GW detectors using the latest models for their cosmological
evolution, and implementing inspiral-merger-ringdown (IMR) gravitational
waveform models in our signal-to-noise ratio calculations. We find that: (1)
the inclusion of the merger/ringdown portion of the signal does not
significantly affect rates for NS-NS and BH-NS systems, but it boosts rates by
a factor $\sim 1.5$ for BH-BH systems; (2) in almost all of our models BH-BH
systems yield by far the largest rates, followed by NS-NS and BH-NS systems,
respectively, and (3) a majority of the detectable BH-BH systems were formed in
the early Universe in low-metallicity environments. We make predictions for the
distributions of detected binaries and discuss what the first GW detections
will teach us about the astrophysics underlying binary formation and evolution.
**URL:** No URL available

## Fast Large Volume Simulations of the 21 cm Signal from the Reionization and pre-Reionization Epochs
**Published:** 2009-11-11
**Authors:** M. B. Silva, M. G. Santos, A. Cooray, L. Ferramacho, A. Amblard
**Abstract:** While limited to low spatial resolution, the next generation low-frequency
radio interferometers that target 21 cm observations during the era of
reionization and prior will have instantaneous fields-of-view that are many
tens of square degrees on the sky. Predictions related to various statistical
measurements of the 21 cm brightness temperature must then be pursued with
numerical simulations of reionization with correspondingly large volume box
sizes, of order 1000 Mpc on one side. We pursue a semi-numerical scheme to
simulate the 21 cm signal during and prior to Reionization by extending a
hybrid approach where simulations are performed by first laying down the linear
dark matter density field, accounting for the non-linear evolution of the
density field based on second-order linear perturbation theory as specified by
the Zel'dovich approximation, and then specifying the location and mass of
collapsed dark matter halos using the excursion-set formalism. The location of
ionizing sources and the time evolving distribution of ionization field is also
specified using an excursion-set algorithm. We account for the brightness
temperature evolution through the coupling between spin and gas temperature due
to collisions, radiative coupling in the presence of Lyman-alpha photons and
heating of the intergalactic medium, such as due to a background of X-ray
photons. The hybrid simulation method we present is capable of producing the
required large volume simulations with adequate resolution in a reasonable time
so a large number of realizations can be obtained with variations in
assumptions related to astrophysics and background cosmology that govern the 21
cm signal.
**URL:** No URL available

## Monte Carlo Method for Calculating Oxygen Abundances and Their Uncertainties from Strong-Line Flux Measurements
**Published:** 2015-05-22
**Authors:** David Fierroz, Yuqian Liu, Federica B. Bianco, Maryam Modjaz, Seung Man Oh, Or Graur, Lisa Kewley
**Abstract:** We present the open-source Python code pyMCZ that determines oxygen abundance
and its distribution from strong emission lines in the standard metallicity
calibrators, based on the original IDL code of Kewley & Dopita (2002) with
updates from Kewley & Ellison (2008), and expanded to include more recently
developed calibrators. The standard strong-line diagnostics have been used to
estimate the oxygen abundance in the interstellar medium through various
emission line ratios in many areas of astrophysics, including galaxy evolution
and supernova host galaxy studies. We introduce a Python implementation of
these methods that, through Monte Carlo sampling, better characterizes the
statistical oxygen abundance confidence region including the effect due to the
propagation of observational uncertainties. These uncertainties are likely to
dominate the error budget in the case of distant galaxies, hosts of cosmic
explosions. Given line flux measurements and their uncertainties, our code
produces synthetic distributions for the oxygen abundance in up to 15
metallicity calibrators simultaneously, as well as for E(B-V), and estimates
their median values and their 68% confidence regions. We test our code on
emission line measurements from a sample of nearby supernova host galaxies (z <
0.15) and compare our metallicity results with those from previous methods. Our
metallicity estimates are consistent with previous methods but yield smaller
statistical uncertainties. Systematic uncertainties are not taken into account.
We offer visualization tools to assess the spread of the oxygen abundance in
the different calibrators, as well as the shape of the estimated oxygen
abundance distribution in each calibrator, and develop robust metrics for
determining the appropriate Monte Carlo sample size. The code is open access
and open source and can be found at https://github.com/nyusngroup/pyMCZ
(Abridged)
**URL:** No URL available

## Introduction to astroML: Machine Learning for Astrophysics
**Published:** 2014-11-18
**Authors:** Jacob T. VanderPlas, Alex Gray, Zeljko Ivezic, Andrew J. Connolly
**Abstract:** Astronomy and astrophysics are witnessing dramatic increases in data volume
as detectors, telescopes and computers become ever more powerful. During the
last decade, sky surveys across the electromagnetic spectrum have collected
hundreds of terabytes of astronomical data for hundreds of millions of sources.
Over the next decade, the data volume will enter the petabyte domain, and
provide accurate measurements for billions of sources. Astronomy and physics
students are not traditionally trained to handle such voluminous and complex
data sets. In this paper we describe astroML; an initiative, based on Python
and scikit-learn, to develop a compendium of machine learning tools designed to
address the statistical needs of the next generation of students and
astronomical surveys. We introduce astroML and present a number of example
applications that are enabled by this package.
**URL:** No URL available

## Machine Learning Interpretability: A Science rather than a tool
**Published:** 2018-07-18
**Authors:** Avinash Mishra, MA Hakim Newton, Abdul Sattar, Abdul Karim
**Abstract:** The term "interpretability" is oftenly used by machine learning researchers
each with their own intuitive understanding of it. There is no universal well
agreed upon definition of interpretability in machine learning. As any type of
science discipline is mainly driven by the set of formulated questions rather
than by different tools in that discipline, e.g. astrophysics is the discipline
that learns the composition of stars, not as the discipline that use the
spectroscopes. Similarly, we propose that machine learning interpretability
should be a discipline that answers specific questions related to
interpretability. These questions can be of statistical, causal and
counterfactual nature. Therefore, there is a need to look into the
interpretability problem of machine learning in the context of questions that
need to be addressed rather than different tools. We discuss about a
hypothetical interpretability framework driven by a question based scientific
approach rather than some specific machine learning model. Using a question
based notion of interpretability, we can step towards understanding the science
of machine learning rather than its engineering. This notion will also help us
understanding any specific problem more in depth rather than relying solely on
machine learning methods.
**URL:** No URL available

## An improved effective-one-body model of spinning, nonprecessing binary black holes for the era of gravitational-wave astrophysics with advanced detectors
**Published:** 2016-11-11
**Authors:** Béla Szilágyi, Heather Fong, Vivien Raymond, Michael Pürrer, Ian Hinder, Ian W. Harry, Lijing Shao, Alejandro Bohé, Tony Chu, Stanislav Babak, Prayush Kumar, Harald P. Pfeiffer, Geoffrey Lovelace, Daniel A. Hemberger, Serguei Ossokine, Andrea Taracchini, Alessandra Buonanno, Mark A. Scheel, Michael Boyle, Lawrence E. Kidder
**Abstract:** We improve the accuracy of the effective-one-body (EOB) waveforms that were
employed during the first observing run of Advanced LIGO for binaries of
spinning, nonprecessing black holes by calibrating them to a set of 141
numerical-relativity (NR) waveforms. The NR simulations expand the domain of
calibration towards larger mass ratios and spins, as compared to the previous
EOBNR model. Merger-ringdown waveforms computed in black-hole perturbation
theory for Kerr spins close to extremal provide additional inputs to the
calibration. For the inspiral-plunge phase, we use a Markov-chain Monte Carlo
algorithm to efficiently explore the calibration space. For the merger-ringdown
phase, we fit the NR signals with phenomenological formulae. After
extrapolation of the calibrated model to arbitrary mass ratios and spins, the
(dominant-mode) EOBNR waveforms have faithfulness --- at design Advanced-LIGO
sensitivity --- above $99\%$ against all the NR waveforms, including 16
additional waveforms used for validation, when maximizing only on initial phase
and time. This implies a negligible loss in event rate due to modeling for
these binary configurations. We find that future NR simulations at mass ratios
$\gtrsim 4$ and double spin $\gtrsim 0.8$ will be crucial to resolve
discrepancies between different ways of extrapolating waveform models. We also
find that some of the NR simulations that already exist in such region of
parameter space are too short to constrain the low-frequency portion of the
models. Finally, we build a reduced-order version of the EOBNR model to speed
up waveform generation by orders of magnitude, thus enabling intensive
data-analysis applications during the upcoming observation runs of Advanced
LIGO.
**URL:** No URL available

## Predictions of the WFIRST Microlensing Survey I: Bound Planet Detection Rates
**Published:** 2018-08-07
**Authors:** Sebastiano Calchi Novati, Nicholas J. Rattenbury, Shude Mao, Matthew T. Penny, Eamonn Kerins, Annie C. Robin, B. Scott Gaudi
**Abstract:** The Wide Field InfraRed Survey Telescope (WFIRST) is the next NASA
astrophysics flagship mission, to follow the James Webb Space Telescope (JWST).
The WFIRST mission was chosen as the top-priority large space mission of the
2010 astronomy and astrophysics decadal survey in order to achieve three
primary goals: to study dark energy via a wide-field imaging survey, to study
exoplanets via a microlensing survey, and to enable a guest observer program.
Here we assess the ability of the several WFIRST designs to achieve the goal of
the microlensing survey to discover a large sample of cold, low-mass exoplanets
with semimajor axes beyond roughly one AU, which are largely impossible to
detect with any other technique. We present the results of a suite of
simulations that span the full range of the proposed WFIRST architectures, from
the original design envisioned by the decadal survey, to the current design,
which utilizes a 2.4-m telescope donated to NASA. By studying such a broad
range of architectures, we are able to determine the impact of design trades on
the expected yields of detected exoplanets. In estimating the yields we take
particular care to ensure that our assumed Galactic model predicts microlensing
event rates that match observations, consider the impact that inaccuracies in
the Galactic model might have on the yields, and ensure that numerical errors
in lightcurve computations do not bias the yields for the smallest mass
exoplanets. For the nominal baseline WFIRST design and a fiducial planet mass
function, we predict that a total of ${\sim}1400$ bound exoplanets with mass
greater than ${\sim}0.1~M_{\oplus}$ should be detected, including ${\sim}200$
with mass ${\lesssim}3~M_{\oplus}$. WFIRST should have sensitivity to planets
with mass down to ${\sim}0.02~M_{\oplus}$, or roughly the mass of Ganymede.
**URL:** No URL available

## Analyzing Inverse Problems with Invertible Neural Networks
**Published:** 2018-08-14
**Authors:** Ullrich Köthe, Carsten Rother, Lena Maier-Hein, Ralf S. Klessen, Daniel Rahner, Sebastian Wirkert, Jakob Kruse, Lynton Ardizzone, Eric W. Pellegrini
**Abstract:** In many tasks, in particular in natural science, the goal is to determine
hidden system parameters from a set of measurements. Often, the forward process
from parameter- to measurement-space is a well-defined function, whereas the
inverse problem is ambiguous: one measurement may map to multiple different
sets of parameters. In this setting, the posterior parameter distribution,
conditioned on an input measurement, has to be determined. We argue that a
particular class of neural networks is well suited for this task -- so-called
Invertible Neural Networks (INNs). Although INNs are not new, they have, so
far, received little attention in literature. While classical neural networks
attempt to solve the ambiguous inverse problem directly, INNs are able to learn
it jointly with the well-defined forward process, using additional latent
output variables to capture the information otherwise lost. Given a specific
measurement and sampled latent variables, the inverse pass of the INN provides
a full distribution over parameter space. We verify experimentally, on
artificial data and real-world problems from astrophysics and medicine, that
INNs are a powerful analysis tool to find multi-modalities in parameter space,
to uncover parameter correlations, and to identify unrecoverable parameters.
**URL:** No URL available

## Multi-timescale analysis of phase transitions in precessing black-hole binaries
**Published:** 2015-06-10
**Authors:** Richard O'Shaughnessy, Ulrich Sperhake, Michael Kesden, Emanuele Berti, Davide Gerosa
**Abstract:** The dynamics of precessing binary black holes (BBHs) in the post-Newtonian
regime has a strong timescale hierarchy: the orbital timescale is very short
compared to the spin-precession timescale which, in turn, is much shorter than
the radiation-reaction timescale on which the orbit is shrinking due to
gravitational-wave emission. We exploit this timescale hierarchy to develop a
multi-scale analysis of BBH dynamics elaborating on the analysis of Kesden et
al. (2015). We solve the spin-precession equations analytically on the
precession time and then implement a quasi-adiabatic approach to evolve these
solutions on the longer radiation-reaction time. This procedure leads to an
innovative "precession-averaged" post-Newtonian approach to studying precessing
BBHs. We use our new solutions to classify BBH spin precession into three
distinct morphologies, then investigate phase transitions between these
morphologies as BBHs inspiral. These precession-averaged post-Newtonian
inspirals can be efficiently calculated from arbitrarily large separations,
thus making progress towards bridging the gap between astrophysics and
numerical relativity.
**URL:** No URL available

## Surprises from the spins: astrophysics and relativity with detections of spinning black-hole mergers
**Published:** 2017-11-27
**Authors:** Davide Gerosa
**Abstract:** Measurements of black-hole spins are of crucial importance to fulfill the
promise of gravitational-wave astronomy. On the astrophysics side, spins are
perhaps the cleanest indicator of black-hole evolutionary processes, thus
providing a preferred way to discriminate how LIGO's black holes form. On the
relativity side, spins are responsible for peculiar dynamical phenomena (from
precessional modulations in the long inspiral to gravitational-wave recoils at
merger) which encode precious information on the underlying astrophysical
processes. I present some examples to explore this deep and fascinating
interplay between spin dynamics (relativity) and environmental effects
(astrophysics). Black-hole spins indeed hide remarkable surprises on both
fronts: morphologies, resonances, constraints on supernova kicks, multiple
merger generations and more...
**URL:** No URL available

## Resonant-plane locking and spin alignment in stellar-mass black-hole binaries: a diagnostic of compact-binary formation
**Published:** 2013-02-18
**Authors:** Richard O'Shaughnessy, Davide Gerosa, Emanuele Berti, Ulrich Sperhake, Michael Kesden
**Abstract:** We study the influence of astrophysical formation scenarios on the
precessional dynamics of spinning black-hole binaries by the time they enter
the observational window of second- and third-generation gravitational-wave
detectors, such as Advanced LIGO/Virgo, LIGO-India, KAGRA and the Einstein
Telescope. Under the plausible assumption that tidal interactions are efficient
at aligning the spins of few-solar mass black-hole progenitors with the orbital
angular momentum, we find that black-hole spins should be expected to
preferentially lie in a plane when they become detectable by gravitational-wave
interferometers. This "resonant plane" is identified by the conditions
\Delta\Phi=0{\deg} or \Delta\Phi=+/-180{\deg}, where \Delta\Phi is the angle
between the components of the black-hole spins in the plane orthogonal to the
orbital angular momentum. If the angles \Delta \Phi can be accurately measured
for a large sample of gravitational-wave detections, their distribution will
constrain models of compact binary formation. In particular, it will tell us
whether tidal interactions are efficient and whether a mechanism such as mass
transfer, stellar winds, or supernovae can induce a mass-ratio reversal (so
that the heavier black hole is produced by the initially lighter stellar
progenitor). Therefore our model offers a concrete observational link between
gravitational-wave measurements and astrophysics. We also hope that it will
stimulate further studies of precessional dynamics, gravitational-wave template
placement and parameter estimation for binaries locked in the resonant plane.
**URL:** No URL available

## ATLAS Probe: Breakthrough Science of Galaxy Evolution, Cosmology, Milky Way, and the Solar System
**Published:** 2018-02-05
**Authors:** Stephen Smee, Daniel Scolnic, Jason Rhodes, Olivier Dore, James Bartlett, Robert Barkhouser, J. Davy Kirkpatrick, Robert Content, Chia-Hsun Chuang, Robert Benjamin, Ranga Chary, Michael J. Hudson, Mario Ballardini, Lauro Moscardini, Jeffrey A. Newman, James Rhoads, Alvaro Orsi, Alice Shapley, Yun Wang, Sangeeta Malhotra, Megan Donahue, Mark Dickinson, Karl Glazebrook, George Helou, Francesco Valentino, Emanuele Daddi, Wesley Fraser, Peter Eisenhardt, Massimo Robberto, Lynne A. Hillenbrand, Charlie Conroy, Zoran Ninkov, Risa H. Wechsler, Peter Behroozi, Michael Ressler, Jarle Brinchmann, Henry C. Ferguson, Christopher Hirata, Andrea Cimatti
**Abstract:** ATLAS (Astrophysics Telescope for Large Area Spectroscopy) Probe is a concept
for a NASA probe-class space mission. It is the follow-up space mission to
WFIRST, boosting its scientific return by obtaining deep IR slit spectroscopy
for 70% of all galaxies imaged by a 2000 sq deg WFIRST High Latitude Survey at
z>0.5. ATLAS will measure accurate and precise redshifts for 200M galaxies out
to z < 7, and deliver spectra that enable a wide range of diagnostic studies of
the physical properties of galaxies over most of cosmic history. ATLAS Probe
science spans four broad categories: (1) Revolutionizing galaxy evolution
studies by tracing the relation between galaxies and dark matter from galaxy
groups to cosmic voids and filaments, from the epoch of reionization through
the peak era of galaxy assembly; (2) Opening a new window into the dark
Universe by weighing the dark matter filaments using 3D weak lensing with
spectroscopic redshifts, and obtaining definitive measurements of dark energy
and modification of General Relativity using galaxy clustering; (3) Probing the
Milky Way's dust-enshrouded regions, reaching the far side of our Galaxy; and
(4) Exploring the formation history of the outer Solar System by characterizing
Kuiper Belt Objects. ATLAS Probe is a 1.5m telescope with a field of view of
0.4 sq deg, and uses Digital Micro-mirror Devices (DMDs) as slit selectors. It
has a spectroscopic resolution of R = 1000 over 1-4 microns, and a
spectroscopic multiplex factor >5,000. ATLAS is designed to fit within the NASA
probe-class space mission cost envelope; it has a single instrument, a
telescope aperture that allows for a lighter launch vehicle, and mature
technology. ATLAS Probe will lead to transformative science over the entire
range of astrophysics: from galaxy evolution to the dark Universe, from Solar
System objects to the dusty regions of the Milky Way.
**URL:** No URL available

## Approximate Bayesian Computation for Forward Modeling in Cosmology
**Published:** 2015-04-27
**Authors:** Alexandre Refregier, Joel Akeret, Sebastian Seehars, Adam Amara, Caspar Hasner
**Abstract:** Bayesian inference is often used in cosmology and astrophysics to derive
constraints on model parameters from observations. This approach relies on the
ability to compute the likelihood of the data given a choice of model
parameters. In many practical situations, the likelihood function may however
be unavailable or intractable due to non-gaussian errors, non-linear
measurements processes, or complex data formats such as catalogs and maps. In
these cases, the simulation of mock data sets can often be made through forward
modeling. We discuss how Approximate Bayesian Computation (ABC) can be used in
these cases to derive an approximation to the posterior constraints using
simulated data sets. This technique relies on the sampling of the parameter
set, a distance metric to quantify the difference between the observation and
the simulations and summary statistics to compress the information in the data.
We first review the principles of ABC and discuss its implementation using a
Population Monte-Carlo (PMC) algorithm and the Mahalanobis distance metric. We
test the performance of the implementation using a Gaussian toy model. We then
apply the ABC technique to the practical case of the calibration of image
simulations for wide field cosmological surveys. We find that the ABC analysis
is able to provide reliable parameter constraints for this problem and is
therefore a promising technique for other applications in cosmology and
astrophysics. Our implementation of the ABC PMC method is made available via a
public code release.
**URL:** No URL available

## astroABC: An Approximate Bayesian Computation Sequential Monte Carlo sampler for cosmological parameter estimation
**Published:** 2016-08-26
**Authors:** Maeve Madigan, Elise Jennings
**Abstract:** Given the complexity of modern cosmological parameter inference where we are
faced with non-Gaussian data and noise, correlated systematics and multi-probe
correlated data sets, the Approximate Bayesian Computation (ABC) method is a
promising alternative to traditional Markov Chain Monte Carlo approaches in the
case where the Likelihood is intractable or unknown. The ABC method is called
"Likelihood free" as it avoids explicit evaluation of the Likelihood by using a
forward model simulation of the data which can include systematics. We
introduce astroABC, an open source ABC Sequential Monte Carlo (SMC) sampler for
parameter estimation. A key challenge in astrophysics is the efficient use of
large multi-probe datasets to constrain high dimensional, possibly correlated
parameter spaces. With this in mind astroABC allows for massive parallelization
using MPI, a framework that handles spawning of jobs across multiple nodes. A
key new feature of astroABC is the ability to create MPI groups with different
communicators, one for the sampler and several others for the forward model
simulation, which speeds up sampling time considerably. For smaller jobs the
Python multiprocessing option is also available. Other key features include: a
Sequential Monte Carlo sampler, a method for iteratively adapting tolerance
levels, local covariance estimate using scikit-learn's KDTree, modules for
specifying optimal covariance matrix for a component-wise or multivariate
normal perturbation kernel, output and restart files are backed up every
iteration, user defined metric and simulation methods, a module for specifying
heterogeneous parameter priors including non-standard prior PDFs, a module for
specifying a constant, linear, log or exponential tolerance level,
well-documented examples and sample scripts. This code is hosted online at
https://github.com/EliseJ/astroABC
**URL:** No URL available

## $γ$-Cascade: A Simple Program to Compute Cosmological Gamma-Ray Propagation
**Published:** 2018-03-30
**Authors:** Carlos Blanco
**Abstract:** Modeling electromagnetic cascades during gamma-ray transport is important in
many applications within astrophysics. This document introduces
{\gamma}-Cascade, a publicly available Mathematica package which allows users
to calculate observed gamma-ray fluxes from point sources as well as from a
distribution of sources. {\gamma}-Cascade semi-analytically computes the
effects of electromagnetic interactions during gamma-ray transport.
**URL:** No URL available

## Learning to Predict the Cosmological Structure Formation
**Published:** 2018-11-15
**Authors:** Barnabás Póczos, Shirley Ho, Siyu He, Siamak Ravanbakhsh, Wei Chen, Yu Feng, Yin Li
**Abstract:** Matter evolved under influence of gravity from minuscule density fluctuations. Non-perturbative structure formed hierarchically over all scales, and developed non-Gaussian features in the Universe, known as the Cosmic Web. To fully understand the structure formation of the Universe is one of the holy grails of modern astrophysics. Astrophysicists survey large volumes of the Universe and employ a large ensemble of computer simulations to compare with the observed data in order to extract the full information of our own Universe. However, to evolve trillions of galaxies over billions of years even with the simplest physics is a daunting task. We build a deep neural network, the Deep Density Displacement Model (hereafter D$^3$M), to predict the non-linear structure formation of the Universe from simple linear perturbation theory. Our extensive analysis, demonstrates that D$^3$M outperforms the second order perturbation theory (hereafter 2LPT), the commonly used fast approximate simulation method, in point-wise comparison, 2-point correlation, and 3-point correlation. We also show that D$^3$M is able to accurately extrapolate far beyond its training data, and predict structure formation for significantly different cosmological parameters. Our study proves, for the first time, that deep learning is a practical and accurate alternative to approximate simulations of the gravitational structure formation of the Universe.
**URL:** No URL available

## Physics, Astrophysics and Cosmology with Gravitational Waves
**Published:** 2009-03-02
**Authors:** B. F. Schutz, B. S. Sathyaprakash
**Abstract:** Gravitational wave detectors are already operating at interesting sensitivity
levels, and they have an upgrade path that should result in secure detections
by 2014. We review the physics of gravitational waves, how they interact with
detectors (bars and interferometers), and how these detectors operate. We study
the most likely sources of gravitational waves and review the data analysis
methods that are used to extract their signals from detector noise. Then we
consider the consequences of gravitational wave detections and observations for
physics, astrophysics, and cosmology.
**URL:** No URL available

## Improving \textsl{Gaia} parallax precision with a data-driven model of stars
**Published:** 2017-06-15
**Authors:** Adrian M. Price-Whelan, David W. Hogg, Lauren Anderson, Jo Bovy, Boris Leistedt
**Abstract:** Converting a noisy parallax measurement into a posterior belief over distance
requires inference with a prior. Usually this prior represents beliefs about
the stellar density distribution of the Milky Way. However, multi-band
photometry exists for a large fraction of the \textsl{\small{Gaia}}
\textsl{\small{TGAS}} Catalog and is incredibly informative about stellar
distances. Here we use \textsl{\small{2MASS}} colors for 1.4 million
\textsl{\small{TGAS}} stars to build a noise-deconvolved empirical prior
distribution for stars in color--magnitude space. This model contains no
knowledge of stellar astrophysics or the Milky Way, but is precise because it
accurately generates a large number of noisy parallax measurements under an
assumption of stationarity; that is, it is capable of combining the information
from many stars. We use the Extreme Deconvolution (\textsl{\small{XD}})
algorithm---an Empirical Bayes approximation to a full hierarchical model of
the true parallax and photometry of every star---to construct this prior. The
prior is combined with a \textsl{\small{TGAS}} likelihood to infer a precise
photometric parallax estimate and uncertainty (and full posterior) for every
star. Our parallax estimates are more precise than the \textsl{\small{TGAS}}
catalog entries by a median factor of 1.2 (14% are more precise by a factor >2)
and are more precise than previous Bayesian distance estimates that use spatial
priors. We validate our parallax inferences using members of the Milky Way star
cluster M67, which is not visible as a cluster in the \textsl{\small{TGAS}}
parallax estimates, but appears as a cluster in our posterior parallax
estimates. Our results, including a parallax posterior pdf for each of 1.4
million \textsl{\small{TGAS}} stars, are available in companion electronic
tables.
**URL:** No URL available

## Exploring galaxy evolution with generative models
**Published:** 2018-12-03
**Authors:** Ce Zhang, Kevin Schawinski, M. Dennis Turp
**Abstract:** Context. Generative models open up the possibility to interrogate scientific
data in a more data-driven way. Aims: We propose a method that uses generative
models to explore hypotheses in astrophysics and other areas. We use a neural
network to show how we can independently manipulate physical attributes by
encoding objects in latent space. Methods: By learning a latent space
representation of the data, we can use this network to forward model and
explore hypotheses in a data-driven way. We train a neural network to generate
artificial data to test hypotheses for the underlying physical processes.
Results: We demonstrate this process using a well-studied process in
astrophysics, the quenching of star formation in galaxies as they move from
low-to high-density environments. This approach can help explore astrophysical
and other phenomena in a way that is different from current methods based on
simulations and observations.
**URL:** No URL available

## Emulation of reionization simulations for Bayesian inference of astrophysics parameters using neural networks
**Published:** 2017-07-31
**Authors:** Jonathan R Pritchard, Claude J Schmit
**Abstract:** Next generation radio experiments such as LOFAR, HERA and SKA are expected to
probe the Epoch of Reionization and claim a first direct detection of the
cosmic 21cm signal within the next decade. Data volumes will be enormous and
can thus potentially revolutionize our understanding of the early Universe and
galaxy formation. However, numerical modelling of the Epoch of Reionization can
be prohibitively expensive for Bayesian parameter inference and how to
optimally extract information from incoming data is currently unclear.
Emulation techniques for fast model evaluations have recently been proposed as
a way to bypass costly simulations. We consider the use of artificial neural
networks as a blind emulation technique. We study the impact of training
duration and training set size on the quality of the network prediction and the
resulting best fit values of a parameter search. A direct comparison is drawn
between our emulation technique and an equivalent analysis using 21CMMC. We
find good predictive capabilities of our network using training sets of as low
as 100 model evaluations, which is within the capabilities of fully numerical
radiative transfer codes.
**URL:** No URL available

## Heuristics for Efficient Sparse Blind Source Separation
**Published:** 2018-12-17
**Authors:** Christophe Kervazo, Cecile Chenot, Jerome Bobin
**Abstract:** Sparse Blind Source Separation (sparse BSS) is a key method to analyze
multichannel data in fields ranging from medical imaging to astrophysics.
However, since it relies on seeking the solution of a non-convex penalized
matrix factorization problem, its performances largely depend on the
optimization strategy. In this context, Proximal Alternating Linearized
Minimization (PALM) has become a standard algorithm which, despite its
theoretical grounding, generally provides poor practical separation results. In
this work, we propose a novel strategy that combines a heuristic approach with
PALM. We show its relevance on realistic astrophysical data.
**URL:** No URL available

## Rank-Approximate Nearest Neighbor Search: Retaining Meaning and Speed in High Dimensions
**Published:** 2009-12-01
**Authors:** Alexander G. Gray, Parikshit Ram, Dongryeol Lee, Hua Ouyang
**Abstract:** The long-standing problem of efficient nearest-neighbor (NN) search has ubiquitous applications ranging from astrophysics to MP3 fingerprinting to bioinformatics to movie recommendations.  As the dimensionality of the dataset increases, exact NN search becomes computationally prohibitive; (1+eps)-distance-approximate NN search can provide large speedups but risks losing the meaning of NN search present in the ranks (ordering) of the distances. This paper presents a simple, practical algorithm allowing the user to, for the first time, directly control the true accuracy of NN search (in terms of ranks) while still achieving the large speedups over exact NN.  Experiments with high-dimensional datasets show that it often achieves faster and more accurate results than the best-known distance-approximate method, with much more stable behavior.
**URL:** No URL available

## A Simple Algorithm for Scalable Monte Carlo Inference
**Published:** 2019-01-02
**Authors:** Alexander Borisenko, Maksym Byshkin, Alessandro Lomi
**Abstract:** The methods of statistical physics are widely used for modelling complex networks. Building on the recently proposed Equilibrium Expectation approach, we derive a simple and efficient algorithm for maximum likelihood estimation (MLE) of parameters of exponential family distributions - a family of statistical models, that includes Ising model, Markov Random Field and Exponential Random Graph models. Computational experiments and analysis of empirical data demonstrate that the algorithm increases by orders of magnitude the size of network data amenable to Monte Carlo based inference. We report results suggesting that the applicability of the algorithm may readily be extended to the analysis of large samples of dependent observations commonly found in biology, sociology, astrophysics, and ecology.
**URL:** No URL available

## Automated Prototype for Asteroids Detection
**Published:** 2019-01-29
**Authors:** O. Vaduvescu, D. Copandean, D. Gorgan
**Abstract:** Near Earth Asteroids (NEAs) are discovered daily, mainly by few major
surveys, nevertheless many of them remain unobserved for years, even decades.
Even so, there is room for new discoveries, including those submitted by
smaller projects and amateur astronomers. Besides the well-known surveys that
have their own automated system of asteroid detection, there are only a few
software solutions designed to help amateurs and mini-surveys in NEAs
discovery. Some of these obtain their results based on the blink method in
which a set of reduced images are shown one after another and the astronomer
has to visually detect real moving objects in a series of images. This
technique becomes harder with the increase in size of the CCD cameras. Aiming
to replace manual detection we propose an automated pipeline prototype for
asteroids detection, written in Python under Linux, which calls some 3rd party
astrophysics libraries.
**URL:** No URL available

## Towards Machine-assisted Meta-Studies: The Hubble Constant
**Published:** 2019-01-31
**Authors:** Rupert A. C. Croft, Thomas D. Kitching, Sebastian Riedel, Daisuke Kawata, Tom Crossland, Pontus Stenetorp
**Abstract:** We present an approach for automatic extraction of measured values from the astrophysical literature, using the Hubble constant for our pilot study. Our rules-based model -- a classical technique in natural language processing -- has successfully extracted 298 measurements of the Hubble constant, with uncertainties, from the 208,541 available arXiv astrophysics papers. We have also created an artificial neural network classifier to identify papers in arXiv which report novel measurements. From the analysis of our results we find that reporting measurements with uncertainties and the correct units is critical information when distinguishing novel measurements in free text. Our results correctly highlight the current tension for measurements of the Hubble constant and recover the $3.5\sigma$ discrepancy -- demonstrating that the tool presented in this paper is useful for meta-studies of astrophysical measurements from a large number of publications.
**URL:** No URL available

## Deep Learning for Multi-Messenger Astrophysics: A Gateway for Discovery in the Big Data Era
**Published:** 2019-02-01
**Authors:** JinJun Xiong, Brigitta M. Sipőcz, Hongyu Shen, Stuart L. Shapiro, Aaron Saxton, Milton Ruiz, Kenton McHenry, Ashish Mahabal, Xin Liu, Roland Haas, Matias Carrasco Kind, G. Bruce Berriman, Zachariah B. Etienne, William Gropp, Timothy J. Williams, Lunan Sun, John Towns, Ed Seidel, Asad Khan, Alexander R. Olivas Jr, Zhizhen Zhao, Philip S. Cowperthwaite, Federica B. Bianco, Elise Jennings, Anushri Gupta, William T. C. Kramer, Wei Wei, Tom Gibbs, Steve Oberlin, M. S. Neubauer, Minsik Cho, Matthew Graham, Igor Andreoni, Gabrielle Allen, E. A. Huerta, Daniel S. Katz, Bernard Schutz, Alex Schwing, Yue Shen, Volodymyr Kindratenko, Shawn Rosofsky, Rahul Biswas, Kyle Chard, J. M. Miller, Jack Wells, Etienne Bachelet, Daniel George, Antonios Tsokaros
**Abstract:** This report provides an overview of recent work that harnesses the Big Data
Revolution and Large Scale Computing to address grand computational challenges
in Multi-Messenger Astrophysics, with a particular emphasis on real-time
discovery campaigns. Acknowledging the transdisciplinary nature of
Multi-Messenger Astrophysics, this document has been prepared by members of the
physics, astronomy, computer science, data science, software and
cyberinfrastructure communities who attended the NSF-, DOE- and NVIDIA-funded
"Deep Learning for Multi-Messenger Astrophysics: Real-time Discovery at Scale"
workshop, hosted at the National Center for Supercomputing Applications,
October 17-19, 2018. Highlights of this report include unanimous agreement that
it is critical to accelerate the development and deployment of novel,
signal-processing algorithms that use the synergy between artificial
intelligence (AI) and high performance computing to maximize the potential for
scientific discovery with Multi-Messenger Astrophysics. We discuss key aspects
to realize this endeavor, namely (i) the design and exploitation of scalable
and computationally efficient AI algorithms for Multi-Messenger Astrophysics;
(ii) cyberinfrastructure requirements to numerically simulate astrophysical
sources, and to process and interpret Multi-Messenger Astrophysics data; (iii)
management of gravitational wave detections and triggers to enable
electromagnetic and astro-particle follow-ups; (iv) a vision to harness future
developments of machine and deep learning and cyberinfrastructure resources to
cope with the scale of discovery in the Big Data Era; (v) and the need to build
a community that brings domain experts together with data scientists on equal
footing to maximize and accelerate discovery in the nascent field of
Multi-Messenger Astrophysics.
**URL:** No URL available

## On the Kelvin-Helmholtz instability with smooth initial conditions -- Linear theory and simulations
**Published:** 2019-02-04
**Authors:** Christoph Pfrommer, Thomas Berlok
**Abstract:** The Kelvin-Helmholtz instability (KHI) is a standard test of hydrodynamic and
magnetohydrodynamic (MHD) simulation codes and finds many applications in
astrophysics. The classic linear theory considers a discontinuity in density
and velocity at the interface of two fluids. However, for numerical simulations
of the KHI such initial conditions do not yield converged results even at the
linear stage of the instability. Instead, smooth profiles of velocity and
density are required for convergence. This renders analytical theory to be only
approximately valid and hinders quantitative comparisons between the classical
theory and simulations. In this paper we derive a linear theory for the KHI
with smooth profiles and illustrate code testing with the MHD code Athena. We
provide the linear solution for the KHI with smooth initial conditions in three
different limits: inviscid hydrodynamics, ideal MHD and Braginskii-MHD. These
linear solutions are obtained numerically with the framework Psecas
(Pseudo-Spectral Eigenvalue Calculator with an Automated Solver), which
generates and solves numerical eigenvalue problems using an equation-parser and
pseudo-spectral methods. The Athena simulations are carried out on a periodic,
Cartesian domain which is useful for code testing purposes. Using Psecas and
analytic theory, we outline the differences between this artificial numerical
setup and the KHI on an infinite Cartesian domain and the KHI in cylindrical
geometry. We discuss several astrophysical applications, such as cold flows in
galaxy formation and cold fronts in galaxy cluster mergers. Psecas, and the
linear solutions used for code testing, are publicly available and can be
downloaded from the web.
**URL:** No URL available

## Extragalactic Searches for Dark Matter Annihilation
**Published:** 2018-09-12
**Authors:** Siddharth Mishra-Sharma
**Abstract:** We are at the dawn of a data-driven era in astrophysics and cosmology. A
large number of ongoing and forthcoming experiments combined with an
increasingly open approach to data availability offer great potential in
unlocking some of the deepest mysteries of the Universe. Among these is
understanding the nature of dark matter (DM)---one of the major unsolved
problems in particle physics. Characterizing DM through its astrophysical
signatures will require a robust understanding of its distribution in the sky
and the use of novel statistical methods.
  The first part of this thesis describes the implementation of a novel
statistical technique which leverages the "clumpiness" of photons originating
from point sources (PSs) to derive the properties of PS populations hidden in
astrophysical datasets. This is applied to data from the Fermi satellite at
high latitudes ($|b| > 30$\deg) to characterize the contribution of PSs of
extragalactic origin. We find that the majority of extragalactic gamma-ray
emission can be ascribed to unresolved PSs having properties consistent with
known sources such as active galactic nuclei. This leaves considerably less
room for significant dark matter contribution.
  The second part of this thesis poses the question: "what is the best way to
look for annihilating dark matter in extragalactic sources?" and attempts to
answer it by constructing a pipeline to robustly map out the distribution of
dark matter outside the Milky Way using galaxy group catalogs. This framework
is then applied to Fermi data and existing group catalogs to search for
annihilating dark matter in extragalactic galaxies and clusters.
**URL:** No URL available

## A Survey of Crowdsourcing in Medical Image Analysis
**Published:** 2019-02-25
**Authors:** Silas Ørting, Panagiotis Mavridis, Helen Spiers, Arno van Hilten, Veronika Cheplygina, Andrew Doyle, Oana Inel, Matthias Hirth, Christopher R. Madan
**Abstract:** Rapid advances in image processing capabilities have been seen across many domains, fostered by the application of machine learning algorithms to "big-data". However, within the realm of medical image analysis, advances have been curtailed, in part, due to the limited availability of large-scale, well-annotated datasets. One of the main reasons for this is the high cost often associated with producing large amounts of high-quality meta-data. Recently, there has been growing interest in the application of crowdsourcing for this purpose; a technique that has proven effective for creating large-scale datasets across a range of disciplines, from computer vision to astrophysics. Despite the growing popularity of this approach, there has not yet been a comprehensive literature review to provide guidance to researchers considering using crowdsourcing methodologies in their own medical imaging analysis. In this survey, we review studies applying crowdsourcing to the analysis of medical images, published prior to July 2018. We identify common approaches, challenges and considerations, providing guidance of utility to researchers adopting this approach. Finally, we discuss future opportunities for development within this emerging domain.
**URL:** No URL available

## JSPAM: A restricted three-body code for simulating interacting galaxies
**Published:** 2015-11-16
**Authors:** Allen Harvey, Anthony Holincheck, John Wallin
**Abstract:** Restricted three-body codes have a proven ability to recreate much of the
disturbed morphology of actual interacting galaxies. As more sophisticated
n-body models were developed and computer speed increased, restricted
three-body codes fell out of favor. However, their supporting role for
performing wide searches of parameter space when fitting orbits to real systems
demonstrates a continuing need for their use. Here we present the model and
algorithm used in the JSPAM code. A precursor of this code was originally
described in 1990, and was called SPAM. We have recently updated the software
with an alternate potential and a treatment of dynamical friction to more
closely mimic the results from n-body tree codes. The code is released publicly
for use under the terms of the Academic Free License (AFL) v.3.0 and has been
added to the Astrophysics Source Code Library.
**URL:** No URL available

## A Parallel Data Compression Framework for Large Scale 3D Scientific Data
**Published:** 2019-03-18
**Authors:** Panagiotis Hadjidoukas, Fabian Wermelinger
**Abstract:** Large scale simulations of complex systems ranging from climate and
astrophysics to crowd dynamics, produce routinely petabytes of data and are
projected to reach the zettabytes level in the coming decade. These simulations
enable unprecedented insights but at the same their effectiveness is hindered
by the enormous data sizes associated with the computational elements and
respective output quantities of interest that impose severe constraints on
storage and I/O time. In this work, we address these challenges through a novel
software framework for scientific data compression. The software (CubismZ)
incorporates efficient wavelet based techniques and the state-of-the-art ZFP,
SZ and FPZIP floating point compressors. The framework relies on a
block-structured data layout, benefits from OpenMP and MPI and targets
supercomputers based on multicores. CubismZ can be used as a tool for ex situ
(offline) compression of scientific datasets and supports conventional
Computational Fluid Dynamics (CFD) file formats. Moreover, it provides a
testbed of comparison, in terms of compression factor and peak signal-to-noise
ratio, for a number of available data compression methods. The software yields
in situ compression ratios of 100x or higher for fluid dynamics data produced
by petascale simulations of cloud cavitation collapse using
$\mathcal{O}(10^{11})$ grid cells, with negligible impact on the total
simulation time.
**URL:** No URL available

## Velocity-induced Acoustic Oscillations at Cosmic Dawn
**Published:** 2019-04-16
**Authors:** Julian B. Muñoz
**Abstract:** The redshifted 21-cm line of hydrogen holds great potential for the study of
cosmology, as it can probe otherwise unobservable cosmic epochs. In particular,
measurements of the 21-cm power spectrum during cosmic dawn---the era when the
first stars were formed---will provide us with a wealth of information about
the astrophysics of stellar formation, as well as the origin of fluctuations in
our Universe. In addition to their usually considered density fluctuations,
dark matter and baryons possess large relative velocities with respect to each
other, due to the baryon acoustic oscillations (BAOs) suffered by the latter,
which suppress the formation of stars during cosmic dawn, leaving an imprint on
21-cm observables during this era. Here we present 21cmvFAST, a version of the
publicly available code 21cmFAST modified to account for this effect. Previous
work has shown that the inclusion of relative velocities produces an acoustic
modulation on the large-scale 21-cm power spectrum during cosmic dawn. By
comparing analytic calculations with simulations from 21cmvFAST, here we
demonstrate that this modulation takes the form of striking velocity-induced
acoustic oscillations (VAOs), both during the Lyman-$\alpha$ coupling era and
the subsequent epoch of heating. The unique shape of these VAOs, which is
determined by the acoustic physics that generated the relative velocities, can
be analytically computed and is easily distinguishable from the usual
density-sourced fluctuations. We find that, under a broad range of
astrophysical assumptions, the VAOs are detectable at high significance by the
upcoming HERA interferometer, which could therefore confirm the presence of
acoustic oscillations at cosmic dawn.
**URL:** No URL available

## A Standard Ruler at Cosmic Dawn
**Published:** 2019-04-16
**Authors:** Julian B. Muñoz
**Abstract:** The matter in our Universe comes in two flavors: dark and baryonic. Of these,
only the latter couples to photons, giving rise to the well-known baryon
acoustic oscillations (BAOs) and, in the process, generating supersonic
relative velocities between dark matter and baryons. These
velocities---imprinted with the acoustic scale in their genesis---impede the
formation of the first stars during cosmic dawn ($z\sim 20$), modulating the
expected 21-cm signal from this era. In a companion paper we showed, combining
numerical simulations and analytic models, that this modulation takes the form
of striking velocity-induced acoustic oscillations (VAOs), with a
well-understood shape that is frozen at recombination, and unaffected by the
unknown astrophysics of star formation. Here we propose using these VAOs as a
standard ruler at cosmic dawn. We find that three years of 21-cm power-spectrum
data from the upcoming HERA interferometer will be able to measure the Hubble
expansion rate $H(z)$ at $z=15-20$ to percent-level precision, ranging from
0.3\% to 11\% depending on the strength of astrophysical feedback processes and
foregrounds. This would provide a new handle on the expansion rate of our
Universe during an otherwise unprobed epoch, opening a window to the mysterious
cosmic-dawn era.
**URL:** No URL available

## Simultaneously constraining the astrophysics of reionisation and the epoch of heating with 21CMMC
**Published:** 2017-05-09
**Authors:** Bradley Greig, Andrei Mesinger
**Abstract:** The cosmic 21 cm signal is set to revolutionise our understanding of the
early Universe, allowing us to probe the 3D temperature and ionisation
structure of the intergalactic medium (IGM). It will open a window onto the
unseen first galaxies, showing us how their UV and X-ray photons drove the
cosmic milestones of the epoch of reionisation (EoR) and epoch of heating
(EoH). To facilitate parameter inference from the 21 cm signal, we previously
developed 21CMMC: a Monte Carlo Markov Chain sampler of 3D EoR simulations.
Here we extend 21CMMC to include simultaneous modelling of the EoH, resulting
in a complete Bayesian inference framework for the astrophysics dominating the
observable epochs of the cosmic 21 cm signal. We demonstrate that second
generation interferometers, the Hydrogen Epoch of Reionisation Array (HERA) and
Square Kilometre Array (SKA) will be able to constrain ionising and X-ray
source properties of the first galaxies with a fractional precision of order
$\sim1$-10 per cent (1$\sigma$). The ionisation history of the Universe can be
constrained to within a few percent. Using our extended framework, we quantify
the bias in EoR parameter recovery incurred by the common simplification of a
saturated spin temperature in the IGM. Depending on the extent of overlap
between the EoR and EoH, the recovered astrophysical parameters can be biased
by $\sim3-10\sigma$.
**URL:** No URL available

## 21CMMC with a 3D light-cone: the impact of the co-evolution approximation on the astrophysics of reionisation and cosmic dawn
**Published:** 2018-01-05
**Authors:** Bradley Greig, Andrei Mesinger
**Abstract:** We extend 21CMMC, a Monte Carlo Markov Chain sampler of 3D reionisation
simulations, to perform parameter estimation directly on 3D light-cones of the
cosmic 21cm signal. This brings theoretical analysis closer to the tomographic
21-cm observations achievable with next generation interferometers like HERA
and the SKA. Parameter recovery can therefore account for modes which evolve
with redshift/frequency. Additionally, simulated data can be more easily
corrupted to resemble real data. Using the light-cone version of 21CMMC, we
quantify the biases in the recovered astrophysical parameters if we use the
21cm power spectrum from the co-evolution approximation to fit a 3D light-cone
mock observation. While ignoring the light-cone effect under most assumptions
will not significantly bias the recovered astrophysical parameters, it can lead
to an underestimation of the associated uncertainty. However significant biases
($\sim$few -- 10 $\sigma$) can occur if the 21cm signal evolves rapidly (i.e.
the epochs of reionisation and heating overlap significantly) and: (i)
foreground removal is very efficient, allowing large physical scales
($k\lesssim0.1$~Mpc$^{-1}$) to be used in the analysis or (ii) theoretical
modelling is accurate to within $\sim10$ per cent in the power spectrum
amplitude.
**URL:** No URL available

## Opportunities in Time-Domain Extragalactic Astrophysics with the NASA Near-Earth Object Camera (NEOCam)
**Published:** 2019-04-12
**Authors:** J. Davy Kirkpatrick, Nicholas P. Ross, Roberto J. Assef, Matthew J. Graham
**Abstract:** This White Pape motivates the time domain extragalactic science case for the
NASA Near-Earth Object Camera (NEOCam). NEOCam is a NASA Planetary mission
whose goal is to discover and characterize asteroids and comets, to assess the
hazard to Earth from near-Earth objects, and to study the origin, evolution,
and fate of asteroids and comets. NEOCam will, however, cover 68% of the
extragalactic sky and as the NEOWISE-R mission has recently proved, infrared
information is now vital for identifying and characterizing the $\gtrsim$10
million IR bright Active Galactic Nuclei, as well as using the IR light curve
to provide deep insights into accretion disk astrophysics. NEOWISE-R data has
also been used to discover Super-luminous Supernovae, dust echos in Tidal
Disruption Events and detects all of the known $z\geq7$ quasars (and over 80%
of the known $z\geq6.70$ quasars). As such, for relatively little additional
cost, adding the capacity for additional NEOCam data processing (and/or
alerting) would have a massive scientific and legacy impact on extragalactic
time domain science.
**URL:** No URL available

## An Ensemble of Bayesian Neural Networks for Exoplanetary Atmospheric Retrieval
**Published:** 2019-05-25
**Authors:** Shawn D. Domagal-Goldman, Atılım Güneş Baydin, Molly D. O'Beirne, Giada N. Arney, Yarin Gal, Frank Soboczenski, Simone Zorzan, Adam D. Cobb, Michael D. Himes, Daniel Angerhausen
**Abstract:** Machine learning is now used in many areas of astrophysics, from detecting exoplanets in Kepler transit signals to removing telescope systematics. Recent work demonstrated the potential of using machine learning algorithms for atmospheric retrieval by implementing a random forest to perform retrievals in seconds that are consistent with the traditional, computationally-expensive nested-sampling retrieval method. We expand upon their approach by presenting a new machine learning model, \texttt{plan-net}, based on an ensemble of Bayesian neural networks that yields more accurate inferences than the random forest for the same data set of synthetic transmission spectra. We demonstrate that an ensemble provides greater accuracy and more robust uncertainties than a single model. In addition to being the first to use Bayesian neural networks for atmospheric retrieval, we also introduce a new loss function for Bayesian neural networks that learns correlations between the model outputs. Importantly, we show that designing machine learning models to explicitly incorporate domain-specific knowledge both improves performance and provides additional insight by inferring the covariance of the retrieved atmospheric parameters. We apply \texttt{plan-net} to the Hubble Space Telescope Wide Field Camera 3 transmission spectrum for WASP-12b and retrieve an isothermal temperature and water abundance consistent with the literature. We highlight that our method is flexible and can be expanded to higher-resolution spectra and a larger number of atmospheric parameters.
**URL:** No URL available

## Generative deep fields: arbitrarily sized, random synthetic astronomical images through deep learning
**Published:** 2019-04-23
**Authors:** James E. Geach, Michael J. Smith
**Abstract:** Generative Adversarial Networks (GANs) are a class of artificial neural
network that can produce realistic, but artificial, images that resemble those
in a training set. In typical GAN architectures these images are small, but a
variant known as Spatial-GANs (SGANs) can generate arbitrarily large images,
provided training images exhibit some level of periodicity. Deep extragalactic
imaging surveys meet this criteria due to the cosmological tenet of isotropy.
Here we train an SGAN to generate images resembling the iconic Hubble Space
Telescope eXtreme Deep Field (XDF). We show that the properties of 'galaxies'
in generated images have a high level of fidelity with galaxies in the real XDF
in terms of abundance, morphology, magnitude distributions and colours. As a
demonstration we have generated a 7.6-billion pixel 'generative deep field'
spanning 1.45 degrees. The technique can be generalised to any appropriate
imaging training set, offering a new purely data-driven approach for producing
realistic mock surveys and synthetic data at scale, in astrophysics and beyond.
**URL:** No URL available

## Bayesian emulator optimisation for cosmology: application to the Lyman-alpha forest
**Published:** 2018-12-11
**Authors:** Andreu Font-Ribera, Hiranya V. Peiris, Keir K. Rogers, Simeon Bird, Licia Verde, Andrew Pontzen
**Abstract:** The Lyman-alpha forest provides strong constraints on both cosmological
parameters and intergalactic medium astrophysics, which are forecast to improve
further with the next generation of surveys including eBOSS and DESI. As is
generic in cosmological inference, extracting this information requires a
likelihood to be computed throughout a high-dimensional parameter space.
Evaluating the likelihood requires a robust and accurate mapping between the
parameters and observables, in this case the 1D flux power spectrum.
Cosmological simulations enable such a mapping, but due to computational time
constraints can only be evaluated at a handful of sample points; "emulators"
are designed to interpolate between these. The problem then reduces to placing
the sample points such that an accurate mapping is obtained while minimising
the number of expensive simulations required. To address this, we introduce an
emulation procedure that employs Bayesian optimisation of the training set for
a Gaussian process interpolation scheme. Starting with a Latin hypercube
sampling (other schemes with good space-filling properties can be used), we
iteratively augment the training set with extra simulations at new parameter
positions which balance the need to reduce interpolation error while focussing
on regions of high likelihood. We show that smaller emulator error from the
Bayesian optimisation propagates to smaller widths on the posterior
distribution. Even with fewer simulations than a Latin hypercube, Bayesian
optimisation shrinks the 95% credible volume by 90% and, e.g., the 1 sigma
error on the amplitude of small-scale primordial fluctuations by 38%. This is
the first demonstration of Bayesian optimisation applied to large-scale
structure emulation, and we anticipate the technique will generalise to many
other probes such as galaxy clustering, weak lensing and 21cm.
**URL:** No URL available

## IGM-Vis: Analyzing Intergalactic and Circumgalactic Medium Absorption Using Quasar Sightlines in a Cosmic Web Context
**Published:** 2018-12-17
**Authors:** Jasmine Otto, Angus G. Forbes, Joseph N. Burchett, David Abramov, J. Xavier Prochaska, Cassia Artanegara
**Abstract:** We introduce IGM-Vis, a novel astrophysics visualization and data analysis
application for investigating galaxies and the gas that surrounds them in
context with their larger scale environment, the Cosmic Web. Environment is an
important factor in the evolution of galaxies from actively forming stars to
quiescent states with little, if any, discernible star formation activity. The
gaseous halos of galaxies (the circumgalactic medium, or CGM) play a critical
role in their evolution, because the gas necessary to fuel star formation and
any gas expelled from widely observed galactic winds must encounter this
interface region between galaxies and the intergalactic medium (IGM). We
present a taxonomy of tasks typically employed in IGM/CGM studies informed by a
survey of astrophysicists at various career levels, and demonstrate how these
tasks are facilitated via the use of our visualization software. Finally, we
evaluate the effectiveness of IGM-Vis through two in-depth use cases that
depict real-world analysis sessions that use IGM/CGM data.
**URL:** No URL available

## An Online Interactive Geometric Database: Including Exact Solutions of Einstein's Field Equations
**Published:** 2001-11-02
**Authors:** Mustapha Ishak, Kayll Lake
**Abstract:** We describe a new interactive database (GRDB) of geometric objects in the general area of differential geometry. Database objects include, but are not restricted to, exact solutions of Einstein's field equations. GRDB is designed for researchers (and teachers) in applied mathematics, physics and related fields. The flexible search environment allows the database to be useful over a wide spectrum of interests, for example, from practical considerations of neutron star models in astrophysics to abstract space-time classification schemes. The database is built using a modular and object-oriented design and uses several Java technologies (e.g. Applets, Servlets, JDBC). These are platform-independent and well adapted for applications developed to run over the World Wide Web. GRDB is accompanied by a virtual calculator (GRTensorJ), a graphical user interface to the computer algebra system GRTensorII used to perform online coordinate, tetrad or basis calculations. The highly interactive nature of GRDB allows for systematic internal self-checking and a minimization of the required internal records.
**URL:** No URL available

## FaVeST: Fast Vector Spherical Harmonic Transforms
**Published:** 2019-07-31
**Authors:** Ming Li, Yu Guang Wang, Quoc T. Le Gia
**Abstract:** Vector spherical harmonics on the unit sphere of $\mathbb{R}^3$ have broad applications in geophysics, quantum mechanics and astrophysics. In the representation of a tangent vector field, one needs to evaluate the expansion and the Fourier coefficients of vector spherical harmonics. In this paper, we develop fast algorithms (FaVeST) for vector spherical harmonic transforms on these evaluations. The forward FaVeST evaluates the Fourier coefficients and has a computational cost proportional to $N\log \sqrt{N}$ for $N$ number of evaluation points. The adjoint FaVeST which evaluates a linear combination of vector spherical harmonics with a degree up to $\sqrt{M}$ for $M$ evaluation points has cost proportional to $M\log\sqrt{M}$. Numerical examples of simulated tangent fields illustrate the accuracy, efficiency and stability of FaVeST.
**URL:** No URL available

## GW170817: Stringent constraints on neutron-star radii from multimessenger observations and nuclear theory
**Published:** 2019-08-27
**Authors:** Stephanie M. Brown, Badri Krishnan, Sumit Kumar, Duncan A. Brown, Collin D. Capano, Sanjay Reddy, Ingo Tews, Soumi De, Ben Margalit
**Abstract:** The properties of neutron stars are determined by the nature of the matter that they contain. These properties can be constrained by measurements of the star's size. We obtain the most stringent constraints on neutron-star radii to date by combining multimessenger observations of the binary neutron-star merger GW170817 with nuclear theory that best accounts for density-dependent uncertainties in the equation of state. We construct equations of state constrained by chiral effective field theory and marginalize over these using the gravitational-wave observations. Combining this with the electromagnetic observations of the merger remnant that imply the presence of a short-lived hyper-massive neutron star, we find that the radius of a $1.4M_\odot$ neutron star is $R_{1.4M_{\odot}} = 11.0^{+0.9}_{-0.6}~{\rm km}$ ($90\%$ credible interval). This constraint has important implications for dense-matter physics and for astrophysics.
**URL:** No URL available

## Constraining the astrophysics and cosmology from 21cm tomography using deep learning with the SKA
**Published:** 2019-07-17
**Authors:** Sultan Hassan, Sambatra Andrianomena, Caitlin Doughty
**Abstract:** Future Square Kilometre Array (SKA) surveys are expected to generate huge datasets of 21cm maps on cosmological scales from the Epoch of Reionization (EoR). We assess the viability of exploiting machine learning techniques, namely, convolutional neural networks (CNN), to simultaneously estimate the astrophysical and cosmological parameters from 21cm maps from semi-numerical simulations. We further convert the simulated 21cm maps into SKA-like mock maps using the detailed SKA antennae distribution, thermal noise and a recipe for foreground cleaning. We successfully design two CNN architectures (VGGNet-like and ResNet-like) that are both efficiently able to extract simultaneously three astrophysical parameters, namely the photon escape fraction (f$_{\rm esc}$), the ionizing emissivity power dependence on halo mass ($C_{\rm ion}$) and the ionizing emissivity redshift evolution index ($D_{\rm ion}$), and three cosmological parameters, namely the matter density parameter ($\Omega_{m}$), the dimensionless Hubble constant ($h$), and the matter fluctuation amplitude ($\sigma_{8}$), from 21cm maps at several redshifts. With the presence of noise from SKA, our designed CNNs are still able to recover these astrophysical and cosmological parameters with great accuracy ($R^{2} > 90\%$), improving to $R^{2} > 95\%$ towards low redshift and low neutral fraction values. Our results show that future 21cm observations can play a key role to break degeneracy between models and tightly constrain the astrophysical and cosmological parameters.
**URL:** No URL available

## Model-independent constraints on dark matter annihilation in dwarf spheroidal galaxies
**Published:** 2018-02-11
**Authors:** Kimberly K. Boddy, Jason Kumar, Pearl Sandick, Danny Marfatia
**Abstract:** We present a general, model-independent formalism for determining bounds on the production of photons in dwarf spheroidal galaxies via dark matter annihilation, applicable to any set of assumptions about dark matter particle physics or astrophysics. As an illustration, we analyze gamma-ray data from the Fermi Large Area Telescope to constrain a variety of nonstandard dark matter models, several of which have not previously been studied in the context of dwarf galaxy searches.
**URL:** No URL available

## The Classical Gravitational N-Body Problem
**Published:** 2005-03-28
**Authors:** Douglas C. Heggie
**Abstract:** Let a number, N, of particles interact classically through Newton's Laws of Motion and Newton's inverse square Law of Gravitation. The resulting equations of motion provide an approximate mathematical model with numerous applications in astrophysics, including the motion of the moon and other bodies in the Solar System (planets, asteroids, comets and meteor particles); stars in stellar systems ranging from binary and other multiple stars to star clusters and galaxies; and the motion of dark matter particles in cosmology. For N=1 and N=2 the equations can be solved analytically. The case N=3 provides one of the richest of all unsolved dynamical problems -- the general three-body problem. For problems dominated by one massive body, as in many planetary problems, approximate methods based on perturbation expansions have been developed. In stellar dynamics, astrophysicists have developed numerous numerical and theoretical approaches to the problem for larger values of N, including treatments based on the Boltzmann equation and the Fokker-Planck equation; such N-body systems can also be modelled as self-gravitating gases, and thermodynamic insights underpin much of our qualitative understanding.
**URL:** No URL available

## End-to-end learning of energy-based representations for irregularly-sampled signals and images
**Published:** 2019-10-01
**Authors:** Lucas. Drumetz, François Rousseau, Ronan Fablet
**Abstract:** For numerous domains, including for instance earth observation, medical imaging, astrophysics,..., available image and signal datasets often involve irregular space-time sampling patterns and large missing data rates. These sampling properties may be critical to apply state-of-the-art learning-based (e.g., auto-encoders, CNNs,...), fully benefit from the available large-scale observations and reach breakthroughs in the reconstruction and identification of processes of interest. In this paper, we address the end-to-end learning of representations of signals, images and image sequences from irregularly-sampled data, i.e. when the training data involved missing data. From an analogy to Bayesian formulation, we consider energy-based representations. Two energy forms are investigated: one derived from auto-encoders and one relating to Gibbs priors. The learning stage of these energy-based representations (or priors) involve a joint interpolation issue, which amounts to solving an energy minimization problem under observation constraints. Using a neural-network-based implementation of the considered energy forms, we can state an end-to-end learning scheme from irregularly-sampled data. We demonstrate the relevance of the proposed representations for different case-studies: namely, multivariate time series, 2D images and image sequences.
**URL:** No URL available

## Revised Radii of Kepler Stars and Planets Using Gaia Data Release 2
**Published:** 2018-05-01
**Authors:** Eric Gaidos, Travis A. Berger, Jennifer L. van Saders, Daniel Huber
**Abstract:** One bottleneck for the exploitation of data from the $Kepler$ mission for stellar astrophysics and exoplanet research has been the lack of precise radii and evolutionary states for most of the observed stars. We report revised radii of 177,911 $Kepler$ stars derived by combining parallaxes from $Gaia$ Data Release 2 with the DR25 $Kepler$ Stellar Properties Catalog. The median radius precision is $\approx$ 8%, a typical improvement by a factor of 4-5 over previous estimates for typical $Kepler$ stars. We find that $\approx$ 67% ($\approx$ 120,000) of all $Kepler$ targets are main-sequence stars, $\approx$ 21% ($\approx$ 37,000) are subgiants, and $\approx$ 12% ($\approx$ 21,000) are red giants, demonstrating that subgiant contamination is less severe than some previous estimates and that Kepler targets are mostly main-sequence stars. Using the revised stellar radii, we recalculate the radii for 2123 confirmed and 1922 candidate exoplanets. We confirm the presence of a gap in the radius distribution of small, close-in planets, but find that the gap is mostly limited to incident fluxes $>$ 200$F_\oplus$ and its location may be at a slightly larger radius (closer to $\approx$ 2$R_\oplus$) when compared to previous results. Further, we find several confirmed exoplanets occupying a previously-described "hot super-Earth desert" at high irradiance, show the relation between gas-giant planet radius and incident flux, and establish a bona-fide sample of eight confirmed planets and 30 planet candidates with $R_{\mathrm{p}}$ $<$ 2$R_\oplus$ in circumstellar "habitable zones" (incident fluxes between 0.25--1.50 $F_\oplus$). The results presented here demonstrate the potential for transformative characterization of stellar and exoplanet populations using $Gaia$ data.
**URL:** No URL available

## Distributed filtered hyperinterpolation for noisy data on the sphere
**Published:** 2019-10-06
**Authors:** Ding-Xuan Zhou, Yu Guang Wang, Shao-Bo Lin
**Abstract:** Problems in astrophysics, space weather research and geophysics usually need to analyze noisy big data on the sphere. This paper develops distributed filtered hyperinterpolation for noisy data on the sphere, which assigns the data fitting task to multiple servers to find a good approximation of the mapping of input and output data. For each server, the approximation is a filtered hyperinterpolation on the sphere by a small proportion of quadrature nodes. The distributed strategy allows parallel computing for data processing and model selection and thus reduces computational cost for each server while preserves the approximation capability compared to the filtered hyperinterpolation. We prove quantitative relation between the approximation capability of distributed filtered hyperinterpolation and the numbers of input data and servers. Numerical examples show the efficiency and accuracy of the proposed method.
**URL:** No URL available

## MADHAT: Model-Agnostic Dark Halo Analysis Tool
**Published:** 2019-10-07
**Authors:** Stephen Hill, Jason Kumar, Pearl Sandick, Barmak Shams Es Haghi, Kimberly K. Boddy
**Abstract:** We present the Model-Agnostic Dark Halo Analysis Tool (MADHAT), a numerical tool which implements a Fermi-LAT data-driven, model-independent analysis of gamma-ray emission from dwarf satellite galaxies and dwarf galaxy candidates due to dark matter annihilation, dark matter decay, or other nonstandard or unknown astrophysics. This tool efficiently provides statistical upper bounds on the number of observed photons in excess of the number expected, based on empirical determinations of foregrounds and backgrounds, using a stacked analysis of any selected set of dwarf targets. It also calculates the resulting bounds on the properties of dark matter under any assumptions the user makes regarding dark sector particle physics or astrophysics. As an application, we determine new bounds on Sommerfeld-enhanced dark matter annihilation in a set of eight dwarfs. MADHAT v1.0 includes 58 dwarfs and dwarf candidate targets, and we discuss future planned developments. MADHAT is available and will be maintained at https://github.com/MADHATdm
**URL:** No URL available

## REBOUND: An open-source multi-purpose N-body code for collisional dynamics
**Published:** 2011-10-21
**Authors:** Shang-Fei Liu, Hanno Rein
**Abstract:** REBOUND is a new multi-purpose N-body code which is freely available under an open-source license. It was designed for collisional dynamics such as planetary rings but can also solve the classical N-body problem. It is highly modular and can be customized easily to work on a wide variety of different problems in astrophysics and beyond. REBOUND comes with three symplectic integrators: leap-frog, the symplectic epicycle integrator (SEI) and a Wisdom-Holman mapping (WH). It supports open, periodic and shearing-sheet boundary conditions. REBOUND can use a Barnes-Hut tree to calculate both self-gravity and collisions. These modules are fully parallelized with MPI as well as OpenMP. The former makes use of a static domain decomposition and a distributed essential tree. Two new collision detection modules based on a plane-sweep algorithm are also implemented. The performance of the plane-sweep algorithm is superior to a tree code for simulations in which one dimension is much longer than the other two and in simulations which are quasi-two dimensional with less than one million particles. In this work, we discuss the different algorithms implemented in REBOUND, the philosophy behind the code's structure as well as implementation specific details of the different modules. We present results of accuracy and scaling tests which show that the code can run efficiently on both desktop machines and large computing clusters.
**URL:** No URL available

## Differentiable Strong Lensing: Uniting Gravity and Neural Nets through Differentiable Probabilistic Programming
**Published:** 2019-10-14
**Authors:** Sydney Otten, Christoph Weniger, Adam Coogan, Paul Hofma, Marco Chianese
**Abstract:** The careful analysis of strongly gravitationally lensed radio and optical images of distant galaxies can in principle reveal DM (sub-)structures with masses several orders of magnitude below the mass of dwarf spheroidal galaxies. However, analyzing these images is a complex task, given the large uncertainties in the source and the lens. Here, we leverage and combine three important computer science developments to approach this challenge from a new perspective. (a) Convolutional deep neural networks, which show extraordinary performance in recognizing and predicting complex, abstract correlation structures in images. (b) Automatic differentiation, which forms the technological backbone for training deep neural networks and increasingly permeates `traditional' physics simulations, thus enabling the application of powerful gradient-based parameter inference techniques. (c) Deep probabilistic programming languages, which not only allow the specification of probabilistic programs and automatize the parameter inference step, but also the direct integration of deep neural networks as model components. In the current work, we demonstrate that it is possible to combine a deconvolutional deep neural network trained on galaxy images as source model with a fully-differentiable and exact implementation of the gravitational lensing physics in a single probabilistic model. This does away with hyperparameter tuning for the source model, enables the simultaneous optimization of nearly one hundred source and lens parameters with gradient-based methods, and allows the use of efficient gradient-based Hamiltonian Monte Carlo posterior sampling techniques. We consider this work as one of the first steps in establishing differentiable probabilistic programming techniques in the particle astrophysics community, which have the potential to significantly accelerate and improve many complex data analysis tasks.
**URL:** No URL available

## Well-balanced finite volume schemes for hydrodynamic equations with general free energy
**Published:** 2018-12-03
**Authors:** Chi-Wang Shu, José A. Carrillo, Serafim Kalliadasis, Sergio P. Perez
**Abstract:** Well balanced and free energy dissipative first- and second-order accurate finite volume schemes are proposed for a general class of hydrodynamic systems with linear and nonlinear damping. The natural Liapunov functional of the system, given by its free energy, allows for a characterization of the stationary states by its variation. An analog property at the discrete level enables us to preserve stationary states at machine precision while keeping the dissipation of the discrete free energy. These schemes allow for analysing accurately the stability properties of stationary states in challeging problems such as: phase transitions in collective behavior, generalized Euler-Poisson systems in chemotaxis and astrophysics, and models in dynamic density functional theories; having done a careful validation in a battery of relevant test cases.
**URL:** No URL available

## An Information Theory Approach on Deciding Spectroscopic Follow Ups
**Published:** 2019-11-06
**Authors:** Javiera Astudillo, Karim Pichara, Pavlos Protopapas, Pablo Huijse
**Abstract:** Classification and characterization of variable phenomena and transient phenomena are critical for astrophysics and cosmology. These objects are commonly studied using photometric time series or spectroscopic data. Given that many ongoing and future surveys are in time-domain and given that adding spectra provide further insights but requires more observational resources, it would be valuable to know which objects should we prioritize to have spectrum in addition to time series. We propose a methodology in a probabilistic setting that determines a-priory which objects are worth taking spectrum to obtain better insights, where we focus 'insight' as the type of the object (classification). Objects for which we query its spectrum are reclassified using their full spectrum information. We first train two classifiers, one that uses photometric data and another that uses photometric and spectroscopic data together. Then for each photometric object we estimate the probability of each possible spectrum outcome. We combine these models in various probabilistic frameworks (strategies) which are used to guide the selection of follow up observations. The best strategy depends on the intended use, whether it is getting more confidence or accuracy. For a given number of candidate objects (127, equal to 5% of the dataset) for taking spectra, we improve 37% class prediction accuracy as opposed to 20% of a non-naive (non-random) best base-line strategy. Our approach provides a general framework for follow-up strategies and can be extended beyond classification and to include other forms of follow-ups beyond spectroscopy.
**URL:** No URL available

## Voxel datacubes for 3D visualization in Blender
**Published:** 2016-11-21
**Authors:** Matías Gárate
**Abstract:** The growth of computational astrophysics and complexity of multidimensional datasets evidences the need for new versatile visualization tools for both analysis and presentation of the data. In this work we show how to use the open source software Blender as a 3D visualization tool to study and visualize numerical simulation results, focusing on astrophysical hydrodynamic experiments. With a datacube as input, the software can generate a volume rendering of the 3D data, show the evolution of a simulation in time, and do a fly-around camera animation to highlight the points of interest. We explain the process to import simulation outputs into Blender using the Voxel Data format, and how to set up a visualization scene in the software interface. This method allows scientists to perform a complementary visual analysis of their data, and display their results in an appealing way, both for outreach and science presentations.
**URL:** No URL available

## Corrfunc: Blazing fast correlation functions with AVX512F SIMD Intrinsics
**Published:** 2019-11-15
**Authors:** Lehman H. Garrison, Manodeep Sinha
**Abstract:** Correlation functions are widely used in extra-galactic astrophysics to extract insights into how galaxies occupy dark matter halos and in cosmology to place stringent constraints on cosmological parameters. A correlation function fundamentally requires computing pair-wise separations between two sets of points and then computing a histogram of the separations. Corrfunc is an existing open-source, high-performance software package for efficiently computing a multitude of correlation functions. In this paper, we will discuss the SIMD AVX512F kernels within Corrfunc, capable of processing 16 floats or 8 doubles at a time. The latest manually implemented Corrfunc AVX512F kernels show a speedup of up to $\sim 4\times$ relative to compiler-generated code for double-precision calculations. The AVX512F kernels show $\sim 1.6\times$ speedup relative to the AVX kernels and compare favorably to a theoretical maximum of $2\times$. In addition, by pruning pairs with too large of a minimum possible separation, we achieve a $\sim 5-10\%$ speedup across all the SIMD kernels. Such speedups highlight the importance of programming explicitly with SIMD vector intrinsics for complex calculations that can not be efficiently vectorized by compilers. Corrfunc is publicly available at https://github.com/manodeep/Corrfunc/.
**URL:** No URL available

## W-boson and trident production in TeV--PeV neutrino observatories
**Published:** 2019-10-23
**Authors:** John F. Beacom, Bei Zhou
**Abstract:** Detecting TeV--PeV cosmic neutrinos provides crucial tests of neutrino physics and astrophysics. The statistics of IceCube and the larger proposed IceCube-Gen2 demand calculations of neutrino-nucleus interactions subdominant to deep-inelastic scattering, which is mediated by weak-boson couplings to nuclei. The largest such interactions are W-boson and trident production, which are mediated instead through photon couplings to nuclei. In a companion paper [1], we make the most comprehensive and precise calculations of those interactions at high energies. In this paper, we study their phenomenological consequences. We find that: (1) These interactions are dominated by the production of on-shell W-bosons, which carry most of the neutrino energy, (2) The cross section on water/iron can be as large as 7.5%/14% that of charged-current deep-inelastic scattering, much larger than the quoted uncertainty on the latter, (3) Attenuation in Earth is increased by as much as 15%, (4) W-boson production on nuclei exceeds that through the Glashow resonance on electrons by a factor of $\simeq$ 20 for the best-fit IceCube spectrum, (5) The primary signals are showers that will significantly affect the detection rate in IceCube-Gen2; a small fraction of events give unique signatures that may be detected sooner.
**URL:** No URL available

## Enabling real-time multi-messenger astrophysics discoveries with deep learning
**Published:** 2019-11-26
**Authors:** JinJun Xiong, John Towns, Brigitta M. Sipőcz, Stuart L. Shapiro, Aaron Saxton, Jonah Miller, Zsuzsa Marka, William T. C. Kramer, Robert Gruendl, Matthew Graham, Francisco Förster, Etienne Bachelet, Zachariah B. Etienne, William Gropp, Timothy J. Williams, Sarah Habib, Matias Carrasco, Lunan Sun, Hongyu Shen, Federica Bianco, Erik Katsavounidis, Ed Seidel, Donald Petravick, Ashish Mahabal, Asad Khan, Zhizhen Zhao, Philip S. Cowperthwaite, Leo Singer, Javier M. Antelis, Elise Jennings, Claudia Moreno, Anushri Gupta, Alexander R. Olivas, Wei Wei, Tom Gibbs, Steve Oberlin, Roland Haas, Minsik Cho, Maya Fishbach, Margaret W. G. Johnson, Kenton McHenry, Igor Andreoni, Gabrielle Allen, E. A. Huerta, Daniel S. Katz, Bruce Berriman, Bernard F. Schutz, Alex Schwing, Yue Shen, Xin Liu, Volodymyr Kindratenko, Shawn Rosofsky, Rahul Biswas, Milton Ruiz, Mark Neubauer, Kyle Chard, Jack Wells, Daniel George, Antonios Tsokaros, Adam Rebei
**Abstract:** Multi-messenger astrophysics is a fast-growing, interdisciplinary field that combines data, which vary in volume and speed of data processing, from many different instruments that probe the Universe using different cosmic messengers: electromagnetic waves, cosmic rays, gravitational waves and neutrinos. In this Expert Recommendation, we review the key challenges of real-time observations of gravitational wave sources and their electromagnetic and astroparticle counterparts, and make a number of recommendations to maximize their potential for scientific discovery. These recommendations refer to the design of scalable and computationally efficient machine learning algorithms; the cyber-infrastructure to numerically simulate astrophysical sources, and to process and interpret multi-messenger astrophysics data; the management of gravitational wave detections to trigger real-time alerts for electromagnetic and astroparticle follow-ups; a vision to harness future developments of machine learning and cyber-infrastructure resources to cope with the big-data requirements; and the need to build a community of experts to realize the goals of multi-messenger astrophysics.
**URL:** No URL available

## Modelling cosmic ray electron physics in cosmological smoothed particle hydrodynamics simulation
**Published:** 2019-12-01
**Authors:** Xiaoli Lian, Dongchao Zheng, Weitian Li, Linfeng Xiao, Jiajun Zhang, Dan Hu, Chenxi Shan, Zhenghao Zhu
**Abstract:** Cosmic ray electron (CRE) acceleration and cooling are important physical processes in astrophysics. We develop an approximative framework to treat CRE physics in the parallel smoothed particle hydrodynamics code Gadget-3. In our methodology, the CRE spectrum of each fluid element is approximated by a single power-law distribution with spatially varying amplitude, upper cut-off, lower cut-off, and spectral index. We consider diffusive shock acceleration to be the source of injection, and oppositely the sinking processes is attributed to synchrotron radiation, inverse Compton scatters, and Coulomb scatters. The adiabatic gains and losses are also included. We show that our formalism produces the energy and pressure with an accuracy of $ > 90\%$ for a free cooling CRE spectrum. Both slope and intensity of the radio emission computed from the CRE population given by our method in cosmological hydro-simulation coincide well with observations, and our results also show that relaxed clusters have lower fluxes. Finally, we investigate several impacts of the CRE processes on the cosmological hydro-simulation, we find that: (1) the pressure of the CRE spectrum is very small and can be ignored in hydro-simulation, (2) the impacts of the CRE processes on the gas phase-space state of hydro-simulation is up to $3\%$, (3) the CRE processes induce a $5\%$ influence on the mass function in the mass range $10^{12} -10^{13} h^{-1} M_{\odot}$, (4) The gas temperature of massive galaxy cluster is influenced by the CRE processes up to $\sim 10\%$.
**URL:** No URL available

## Building high accuracy emulators for scientific simulations with deep neural architecture search
**Published:** 2020-01-17
**Authors:** S. M. Vinko, J. Topp-Mugglestone, M. Jarvis, S. Oliver, D. Watson-Parris, M. F. Kasim, L. Deaconu, E. Viezzer, S. Khatiwala, P. Hatfield, G. Gregori, D. H. Froula, J. Korenaga
**Abstract:** Computer simulations are invaluable tools for scientific discovery. However, accurate simulations are often slow to execute, which limits their applicability to extensive parameter exploration, large-scale data analysis, and uncertainty quantification. A promising route to accelerate simulations by building fast emulators with machine learning requires large training datasets, which can be prohibitively expensive to obtain with slow simulations. Here we present a method based on neural architecture search to build accurate emulators even with a limited number of training data. The method successfully accelerates simulations by up to 2 billion times in 10 scientific cases including astrophysics, climate science, biogeochemistry, high energy density physics, fusion energy, and seismology, using the same super-architecture, algorithm, and hyperparameters. Our approach also inherently provides emulator uncertainty estimation, adding further confidence in their use. We anticipate this work will accelerate research involving expensive simulations, allow more extensive parameters exploration, and enable new, previously unfeasible computational discovery.
**URL:** No URL available

## Multi-class Gaussian Process Classification with Noisy Inputs
**Published:** 2020-01-28
**Authors:** Daniel Hernández-Lobato, Eduardo C. Garrido-Merchán, Carlos Villacampa-Calvo, Bryan Zaldivar
**Abstract:** It is a common practice in the machine learning community to assume that the observed data are noise-free in the input attributes. Nevertheless, scenarios with input noise are common in real problems, as measurements are never perfectly accurate. If this input noise is not taken into account, a supervised machine learning method is expected to perform sub-optimally. In this paper, we focus on multi-class classification problems and use Gaussian processes (GPs) as the underlying classifier. Motivated by a data set coming from the astrophysics domain, we hypothesize that the observed data may contain noise in the inputs. Therefore, we devise several multi-class GP classifiers that can account for input noise. Such classifiers can be efficiently trained using variational inference to approximate the posterior distribution of the latent variables of the model. Moreover, in some situations, the amount of noise can be known before-hand. If this is the case, it can be readily introduced in the proposed methods. This prior information is expected to lead to better performance results. We have evaluated the proposed methods by carrying out several experiments, involving synthetic and real data. These include several data sets from the UCI repository, the MNIST data set and a data set coming from astrophysics. The results obtained show that, although the classification error is similar across methods, the predictive distribution of the proposed methods is better, in terms of the test log-likelihood, than the predictive distribution of a classifier based on GPs that ignores input noise.
**URL:** No URL available

## Ensemble Slice Sampling: Parallel, black-box and gradient-free inference for correlated & multimodal distributions
**Published:** 2020-02-14
**Authors:** Florian Beutler, Minas Karamanis
**Abstract:** Slice Sampling has emerged as a powerful Markov Chain Monte Carlo algorithm that adapts to the characteristics of the target distribution with minimal hand-tuning. However, Slice Sampling's performance is highly sensitive to the user-specified initial length scale hyperparameter and the method generally struggles with poorly scaled or strongly correlated distributions. This paper introduces Ensemble Slice Sampling (ESS), a new class of algorithms that bypasses such difficulties by adaptively tuning the initial length scale and utilising an ensemble of parallel walkers in order to efficiently handle strong correlations between parameters. These affine-invariant algorithms are trivial to construct, require no hand-tuning, and can easily be implemented in parallel computing environments. Empirical tests show that Ensemble Slice Sampling can improve efficiency by more than an order of magnitude compared to conventional MCMC methods on a broad range of highly correlated target distributions. In cases of strongly multimodal target distributions, Ensemble Slice Sampling can sample efficiently even in high dimensions. We argue that the parallel, black-box and gradient-free nature of the method renders it ideal for use in scientific fields such as physics, astrophysics and cosmology which are dominated by a wide variety of computationally expensive and non-differentiable models.
**URL:** No URL available

## A computational theoretical approach for mining data on transient events from databases of high energy astrophysics experiments
**Published:** 2020-04-08
**Authors:** Francesco Lazzarotto, Maria Teresa Pazienza, Marco Feroci
**Abstract:** Data on transient events, like GRBs, are often contained in large databases of unstructured data from space experiments, merged with potentially large amount of background or simply undesired information. We present a computational formal model to apply techniques of modern computer science -such as Data Mining (DM) and Knowledge Discovering in Databases (KDD)- to a generic, large database derived from a high energy astrophysics experiment. This method is aimed to search, identify and extract expected information, and maybe to discover unexpected information .
**URL:** No URL available

## Compressed Convolutional LSTM: An Efficient Deep Learning framework to Model High Fidelity 3D Turbulence
**Published:** 2019-02-28
**Authors:** Michael Chertkov, Don Daniel, Daniel Livescu, Arvind Mohan
**Abstract:** High-fidelity modeling of turbulent flows is one of the major challenges in computational physics, with diverse applications in engineering, earth sciences and astrophysics, among many others. The rising popularity of high-fidelity computational fluid dynamics (CFD) techniques like direct numerical simulation (DNS) and large eddy simulation (LES) have made significant inroads into the problem. However, they remain out of reach for many practical three-dimensional flows characterized by extremely large domains and transient phenomena. Therefore designing efficient and accurate data-driven generative approaches to model turbulence is a necessity. We propose a novel training approach for dimensionality reduction and spatio-temporal modeling of the three-dimensional dynamics of turbulence using a combination of Convolutional autoencoder and the Convolutional LSTM neural networks. The quality of the emulated turbulent fields is assessed with rigorous physics-based statistical tests, instead of visual assessments. The results show significant promise in the training methodology to generate physically consistent turbulent flows at a small fraction of the computing resources required for DNS.
**URL:** No URL available

## BART-based inference for Poisson processes
**Published:** 2020-05-16
**Authors:** Axel Gandy, Sarah Filippi, Seth Flaxman, Mauricio Barahona, Emma McCoy, Stamatina Lamprinakou
**Abstract:** The effectiveness of Bayesian Additive Regression Trees (BART) has been demonstrated in a variety of contexts including non-parametric regression and classification. A BART scheme for estimating the intensity of inhomogeneous Poisson processes is introduced. Poisson intensity estimation is a vital task in various applications including medical imaging, astrophysics and network traffic analysis. The new approach enables full posterior inference of the intensity in a non-parametric regression setting. The performance of the novel scheme is demonstrated through simulation studies on synthetic and real datasets up to five dimensions, and the new scheme is compared with alternative approaches.
**URL:** No URL available

## Distributed Learning via Filtered Hyperinterpolation on Manifolds
**Published:** 2020-07-18
**Authors:** Guido Montúfar, Yu Guang Wang
**Abstract:** Learning mappings of data on manifolds is an important topic in contemporary machine learning, with applications in astrophysics, geophysics, statistical physics, medical diagnosis, biochemistry, 3D object analysis. This paper studies the problem of learning real-valued functions on manifolds through filtered hyperinterpolation of input-output data pairs where the inputs may be sampled deterministically or at random and the outputs may be clean or noisy. Motivated by the problem of handling large data sets, it presents a parallel data processing approach which distributes the data-fitting task among multiple servers and synthesizes the fitted sub-models into a global estimator. We prove quantitative relations between the approximation quality of the learned function over the entire manifold, the type of target function, the number of servers, and the number and type of available samples. We obtain the approximation rates of convergence for distributed and non-distributed approaches. For the non-distributed case, the approximation order is optimal.
**URL:** No URL available

## LRP2020: Machine Learning Advantages in Canadian Astrophysics
**Published:** 2019-10-02
**Authors:** L. Perreault-Levasseur, S. Fabbro, Y. Hezaveh, R. Hlozek, K. M. Yi, K. A. Venn, J. Bovy, S. Ravanbakhsh, S. Ellison, L. Spencer, J. Woo, H. Teimoorinia, G. Eadie, JJ. Kavelaars, A Liu
**Abstract:** The application of machine learning (ML) methods to the analysis of astrophysical datasets is on the rise, particularly as the computing power and complex algorithms become more powerful and accessible. As the field of ML enjoys a continuous stream of breakthroughs, its applications demonstrate the great potential of ML, ranging from achieving tens of millions of times increase in analysis speed (e.g., modeling of gravitational lenses or analysing spectroscopic surveys) to solutions of previously unsolved problems (e.g., foreground subtraction or efficient telescope operations). The number of astronomical publications that include ML has been steadily increasing since 2010. With the advent of extremely large datasets from a new generation of surveys in the 2020s, ML methods will become an indispensable tool in astrophysics. Canada is an unambiguous world leader in the development of the field of machine learning, attracting large investments and skilled researchers to its prestigious AI Research Institutions. This provides a unique opportunity for Canada to also be a world leader in the application of machine learning in the field of astrophysics, and foster the training of a new generation of highly skilled researchers.
**URL:** No URL available

## Towards the Development of Entropy-Based Anomaly Detection in an Astrophysics Simulation
**Published:** 2020-09-05
**Authors:** M. Todd Young, Michael Matheson, Drew Schmidt, Bronson Messer
**Abstract:** The use of AI and ML for scientific applications is currently a very exciting and dynamic field. Much of this excitement for HPC has focused on ML applications whose analysis and classification generate very large numbers of flops. Others seek to replace scientific simulations with data-driven surrogate models. But another important use case lies in the combination application of ML to improve simulation accuracy. To that end, we present an anomaly problem which arises from a core-collapse supernovae simulation. We discuss strategies and early successes in applying anomaly detection techniques from machine learning to this scientific simulation, as well as current challenges and future possibilities.
**URL:** No URL available

## Removing Astrophysics in 21 cm maps with Neural Networks
**Published:** 2020-06-25
**Authors:** Francisco Villaescusa-Navarro, Pablo Villanueva-Domingo
**Abstract:** Measuring temperature fluctuations in the 21 cm signal from the Epoch of Reionization and the Cosmic Dawn is one of the most promising ways to study the Universe at high redshifts. Unfortunately, the 21 cm signal is affected by both cosmology and astrophysics processes in a non-trivial manner. We run a suite of 1,000 numerical simulations with different values of the main astrophysical parameters. From these simulations we produce tens of thousands of 21 cm maps at redshifts $10\leq z\leq 20$. We train a convolutional neural network to remove the effects of astrophysics from the 21 cm maps, and output maps of the underlying matter field. We show that our model is able to generate 2D matter fields that not only resemble the true ones visually, but whose statistical properties agree with the true ones within a few percent down to pretty small scales. We demonstrate that our neural network retains astrophysical information, that can be used to constrain the value of the astrophysical parameters. Finally, we use saliency maps to try to understand which features of the 21 cm maps is the network using in order to determine the value of the astrophysical parameters.
**URL:** No URL available

## Measuring the binary black hole mass spectrum with an astrophysically motivated parameterization
**Published:** 2018-01-08
**Authors:** Colm Talbot, Eric Thrane
**Abstract:** Gravitational-wave detections have revealed a previously unknown population of stellar mass black holes with masses above $20\, M_{\odot}$. These observations provide a new way to test models of stellar evolution for massive stars. By considering the astrophysical processes likely to determine the shape of the binary black hole mass spectrum, we construct a parameterized model to capture key spectral features that relate gravitational-wave data to theoretical stellar astrophysics. In particular, we model the signature of pulsational pair-instability supernovae, which are expected to cause all stars with initial mass $100\, M_{\odot}\lesssim M \lesssim 150\, M_{\odot}$ to form $\sim 40\, M_{\odot}$ black holes. This would cause a cut-off in the black hole mass spectrum along with an excess of black holes near $40\, M_{\odot}$. We carry out a simulated data study to illustrate some of the stellar physics that can be inferred using gravitational-wave measurements of binary black holes and demonstrate several such inferences that might be made in the near future. First, we measure the minimum and maximum stellar black hole mass. Second, we infer the presence of a peak due to pair-instability supernovae. Third, we measure the black hole mass ratio distribution. Finally, we show how inadequate models of the black hole mass spectrum lead to biased estimates of the merger rate and the amplitude of the stochastic gravitational-wave background.
**URL:** No URL available

## A Precessing Numerical Relativity Waveform Surrogate Model for Binary Black Holes: A Gaussian Process Regression Approach
**Published:** 2019-03-21
**Authors:** Jonathan Gair, Daniel Williams, James A Clark, Ik Siong Heng, Bhavesh Khamesra
**Abstract:** Gravitational wave astrophysics relies heavily on the use of matched filtering both to detect signals in noisy data from detectors, and to perform parameter estimation on those signals. Matched filtering relies upon prior knowledge of the signals expected to be produced by a range of astrophysical systems, such as binary black holes. These waveform signals can be computed using numerical relativity techniques, where the Einstein field equations are solved numerically, and the signal is extracted from the simulation. Numerical relativity simulations are, however, computationally expensive, leading to the need for a surrogate model which can predict waveform signals in regions of the physical parameter space which have not been probed directly by simulation. We present a method for producing such a surrogate using Gaussian process regression which is trained directly on waveforms generated by numerical relativity. This model returns not just a single interpolated value for the waveform at a new point, but a full posterior probability distribution on the predicted value. This model is therefore an ideal component in a Bayesian analysis framework, through which the uncertainty in the interpolation can be taken into account when performing parameter estimation of signals.
**URL:** No URL available

## AxioNyx: Simulating Mixed Fuzzy and Cold Dark Matter
**Published:** 2020-07-16
**Authors:** Richard Easther, Jens C. Niemeyer, Christoph Behrens, Mateja Gosenca, Bodo Schwabe
**Abstract:** The distinctive effects of fuzzy dark matter are most visible at non-linear galactic scales. We present the first simulations of mixed fuzzy and cold dark matter, obtained with an extended version of the Nyx code. Fuzzy (or ultralight, or axion-like) dark matter dynamics are governed by the comoving Schr\"odinger-Poisson equation. This is evolved with a pseudospectral algorithm on the root grid, and with finite differencing at up to six levels of adaptive refinement. Cold dark matter is evolved with the existing N-body implementation in Nyx. We present the first investigations of spherical collapse in mixed dark matter models, focusing on radial density profiles, velocity spectra and soliton formation in collapsed halos. We find that the effective granule masses decrease in proportion to the fraction of fuzzy dark matter which quadratically suppresses soliton growth, and that a central soliton only forms if the fuzzy dark matter fraction is greater than 10\%. The Nyx framework supports baryonic physics and key astrophysical processes such as star formation. Consequently, AxioNyx will enable increasingly realistic studies of fuzzy dark matter astrophysics.
**URL:** No URL available

## Evolution of the 21 cm signal throughout cosmic history
**Published:** 2008-02-15
**Authors:** Abraham Loeb, Jonathan R. Pritchard
**Abstract:** The potential use of the redshifted 21 cm line from neutral hydrogen for probing the epoch of reionization is motivating the construction of several low-frequency interferometers. There is also much interest in the possibility of constraining the initial conditions from inflation and the nature of the dark matter and dark energy by probing the power-spectrum of density perturbations in three dimensions and on smaller scales than probed by the microwave background anisotropies. Theoretical understanding of the 21 cm signal has been fragmented into different regimes of physical interest. In this paper, we make the first attempt to describe the full redshift evolution of the 21 cm signal between 0<z<300. We include contributions to the 21 cm signal from fluctuations in the gas density, temperature and neutral fraction, as well as the Lyman alpha flux, and allow for a post-reionization signal from damped Ly alpha systems. Our comprehensive analysis provides a useful foundation for optimizing the design of future arrays whose goal is to separate the particle physics from the astrophysics, either by probing the peculiar velocity distortion of the 21 cm power spectrum, or by extending the 21 cm horizon to z > 25 before the first galaxies had formed, or to z < 6 when the residual pockets of hydrogen trace large scale structure.
**URL:** No URL available

## Dynamical computation of the density of states and Bayes factors using nonequilibrium importance sampling
**Published:** 2018-09-28
**Authors:** Eric Vanden-Eijnden, Grant M. Rotskoff
**Abstract:** Nonequilibrium sampling is potentially much more versatile than its equilibrium counterpart, but it comes with challenges because the invariant distribution is not typically known when the dynamics breaks detailed balance. Here, we derive a generic importance sampling technique that leverages the statistical power of configurations transported by nonequilibrium trajectories, and can be used to compute averages with respect to arbitrary target distributions. As a dissipative reweighting scheme, the method can be viewed in relation to the annealed importance sampling (AIS) method and the related Jarzynski equality. Unlike AIS, our approach gives an unbiased estimator, with provably lower variance than directly estimating the average of an observable. We also establish a direct relation between a dynamical quantity, the dissipation, and the volume of phase space, from which we can compute quantities such as the density of states and Bayes factors. We illustrate the properties of estimators relying on this sampling technique in the context of density of state calculations, showing that it scales favorable with dimensionality -- in particular, we show that it can be used to compute the phase diagram of the mean-field Ising model from a single nonequilibrium trajectory. We also demonstrate the robustness and efficiency of the approach with an application to a Bayesian model comparison problem of the type encountered in astrophysics and machine learning.
**URL:** No URL available

## Monte-Carlo simulations of photohadronic processes in astrophysics
**Published:** 1999-03-31
**Authors:** Todor Stanev, R. J. Protheroe, J. P. Rachen, Ralph Engel, A. Muecke
**Abstract:** A new Monte Carlo program for photohadronic interactions of relativistic nucleons with an ambient photon radiation field is presented. The event generator is designed to fulfil typical astrophysical requirements, but can also be used for radiation and background studies at high energy colliders such as LEP2 and HERA, as well as for simulations of photon induced air showers. We consider the full photopion production cross section from the pion production threshold up to high energies. It includes resonance excitation and decay, direct single pion production and diffractive and non-diffractive multiparticle production. The cross section of each individual process is calculated by fitting experimental data, while the kinematics is determined by the underlying particle production process. We demonstrate that our model is capable of reproducing known accelerator data over a wide energy range.
**URL:** No URL available

