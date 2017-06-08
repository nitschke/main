%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coupled class for LaTeX files                 v1.1 Feb-23-2000 %%
%%                                                                %%
%% The regular article class commands still exist and may be used.%%
%%                                                                %%
%% \maketitle make the first page automatically                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions\relax

\LoadClass[twoside,12pt,a4paper]{article}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{times}
\usepackage{url}

\usepackage[plainpages=false]{hyperref}
\hypersetup{%	bookmarks=true,
		bookmarksopen=true,
		bookmarksnumbered=true,
		pdfborder={0 0 0},
		colorlinks=true,
		linkcolor=blue,
		urlcolor=blue	
	}

\RequirePackage{calc}
\RequirePackage{indentfirst}

%% Title
\renewcommand{\maketitle}{%
 \newpage
 \pagenumbering{gobble}%
 \thispagestyle{empty}
 \begin{minipage}[t][1cm][t]{16cm}
  \begin{flushright}
   \fontsize{10}{12}\selectfont
\vspace{-50pt}
{\bf X-DMS 2017}\\
eXtended Discretization MethodS\\
    June 19-21, 2017, Ume\aa{}, Sweden
  \end{flushright}
 \end{minipage}

 \begin{center}%
  \vspace{-30pt}%
  {%3
   \fontsize{14}{16}\selectfont
   \bfseries
    {\@title}
   \par
  }%3
  {%4
\fontsize{12}{12}\selectfont
   \bfseries\@author
   \par
  }%4
  \vspace{-10pt}%
  {
\fontsize{11}{12}\selectfont
   \normalfont\@address
   \par
  }%5
 \end{center}
 \noindent{\it Key words:\space}{\@keywords}
 \par
 \noindent{\bf ABSTRACT}
 \par
}% End of \@maketitle

\def\author{\@ifnextchar[{\@getshortauthor}{\@getauthor}}
\def\@getauthor#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\def\@getshortauthor[#1]#2{\gdef\@shortauthor{#1}\gdef\@author{#2}}

\newcommand{\address}[1]{\gdef\@address{#1}}

\newcommand{\keywords}[1]{\gdef\@keywords{#1}}

\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}

\newcommand{\heading}[1]{\gdef\@heading{#1}}

\def\ps@myheadings{%
    \def\@oddfoot{\normalfont\hfil\thepage\hfil}
    \let\@evenfoot\@oddfoot
      \def\@evenhead{%
        %\parbox[t][11.5pt][t]{\textwidth}{%
\parbox[t][12.5pt][t]{\textwidth}{%
 		\vspace{-10pt}%        
      \centering
       \normalfont
       \fontsize{10}{10}\selectfont
       \@heading\\
\vspace{-.2cm}
       \rule{\textwidth}{0.5pt}%
    }
   }%
    \let\@oddhead\@evenhead
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
}

\def\sectionLevel{1}
\def\sectionIndent{0pt}
\def\sectionSpaceBefore{12pt}
\def\sectionSpaceAfter{6pt}
\def\sectionStyle{\normalsize\bfseries\MakeUppercase}

\renewcommand{\section}{%
  \@startsection{section}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%


\renewcommand{\@makecaption}[2]{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\normalsize\fontsize{12}{12}\selectfont {\bf #1}: #2}
  \ifdim \wd\@tempboxa >\hsize
    \normalsize\fontsize{10}{12}\selectfont #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

\renewenvironment{thebibliography}[1]
  {
  \vspace{-10pt}
   \section*{\MakeUppercase\refname\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
    \vspace{-15pt}

   \list{\@biblabel{\@arabic\c@enumiv}}%
        {\settowidth\labelwidth{\@biblabel{#1}}%
         \leftmargin\labelwidth
         \advance\leftmargin\labelsep
         \@openbib@code
         \usecounter{enumiv}%
         \let\p@enumiv\@empty
         \renewcommand\theenumiv{\@arabic\c@enumiv}}%
   \sloppy
   \clubpenalty4000
   \@clubpenalty \clubpenalty
   \widowpenalty4000%
   \sfcode`\.\@m}
  {\def\@noitemerr
   {\@latex@warning{Empty `thebibliography' environment}}%
   \endlist}
\setlength{\topmargin}{-0.6cm}
\setlength{\evensidemargin}{-5mm}
\setlength{\oddsidemargin}{5mm}
\setlength{\textheight}{230mm}
\setlength{\textwidth}{160mm}
\setlength{\textfloatsep}{12pt}
\setlength{\abovedisplayskip}{6pt}
\setlength{\belowdisplayskip}{6pt}
\setlength{\parindent}{0mm}
\setlength{\parskip}{12pt}
\setlength{\abovecaptionskip}{6pt}
\setlength{\belowcaptionskip}{6pt}
\setlength{\itemsep}{0pt}

\AtBeginDocument{%
  \pagestyle{myheadings}
  \maketitle
  \let\maketitle\relax
}

\endinput
