%% The Evolution of the Genetic Code
%

% This website is a good starting point.
web( 'http://www.evolvingcode.net' );

%% How is DNA translated?
%
% There are many steps required to translate DNA into a protein.  The first
% is to transcribe the target section of DNA into a strand of messenger RNA
% (mRNA).  Without this step, all protein creation would need to take place
% around the DNA, which would result in over-crowding and a severe loss of
% efficiency in the cell.  The ability to create multiple copies of
% particular sections of DNA also allows for control over the quantities of
% proteins produced, depending upon the current requirements.
%
% The next step is to translate the mRNA into a strand of amino acids.
% Translating a sequence of codons (3 nucelotides) into a sequence of amino 
% acids requires some sort of language, otherwise the sequences would make
% no sense.  This language is defined by strands of transfer RNA (tRNA), which are
% themselves a part of the DNA.  tRNA is a non-coding RNA that functions without
% translation into a protein.
%
% Transfer RNA is a short chain, usually less than 100 bases long.  It
% forms a 'clover leaf' shape due to the bonding between short streches of
% complementary bases.  This then folds to form an 'L' shaped structure,
% with the _anticodon_ at one end and an _acceptor stem_ at the other.  The
% acceptor stem binds to the amino acid that is to be coded and the
% anticodon is the reverse complement of the codon that is matches.  For
% example, the codon 'GCA' is matched by anticodon 'TGC' and codes for the
% amino acid Cysteine.
% 
% There is a lot more detail to this process that help to explain issues
% such as how only 20 tRNAs can match 64 codons.  This will be covered
% in later sections.
%
figTrna = imread( 'img/trna.png' );
fig = figure; image( figTrna ); 
axis off;
title( 'A tRNA Schematic (wikipedia.org)' );
annotation( fig,'textarrow',[0.35 0.4853],[0.2404 0.2659],...
    'String',{'Anticodon'} );
annotation( fig,'textarrow',[0.4265 0.5397],[0.8342 0.8051],...
    'String',{'Acceptor Stem'} );

figTrna3d = imread( 'img/3d_TRNA.png' );
fig = figure; image( figTrna3d );
axis off;
title( {'3D Structure of tRNA (wikipedia.org)'} );
annotation( fig,'textarrow',[0.6643 0.6179],[0.5714 0.6857],...
    'String',{'Acceptor Stem'} );
annotation( fig,'textarrow',[0.5375 0.3929],[0.1952 0.2024],...
    'String',{'Anticodon'} );


%% The Genetic Codes
%
% Although the standard genetic code is often referred to as 'universal', there are
% in fact a number of variations.  The bacterial and plant plastid code is
% actually the same as the standard code, except that it uses four extra
% start codons.  There is greater difference with the Yeast Mitochondrial
% code, which differs at six codons, as well as in the starts.
%
% The following function displays the differences between the specified
% pair of genetic codes, in this case the previous examples are shown.
%
codingDifferences( 'Standard', 'Bacterial' );
codingDifferences( 'Standard', 'Yeast Mitochondrial' );

%% A Tree of Genetic Codes
%
% A simple measure of the similarity between genetic codes is one which
% counts the number of times that codes disagree about codon to amino acid
% mapping.
%
[gcDist gcNames] = gcDistance( 'disagreementCount' );
gcTree = seqlinkage( gcDist, 'UPGMA', gcNames );

roots = { 'Standard', 'Trematode' };
for rootID = 1 : length( roots )
    rootName = roots( rootID );
    root = getbyname( gcTree, rootName );
    gcTree = reroot( gcTree, root );
    plot( gcTree, 'type', 'square' );
    title( rootName );
end


%% The Redundancy of the Genetic Code
%
% After the structure of DNA was discovered by James D. Watson and Francis 
% Crick in 1953, there was a great deal of work towards discovering how
% strings made up of four nucleotides could encode 20 amino acids.  

% Initial theories were highly ordered and geared toward efficient storage 
% and transmission [Hayes 1998].  The problem of reading frames was
% understood early on, leading to the idea of 'comma-free' codes.  These
% assign only a single codon to each amino acid, with the remaining 44
% unused codons having no meaning.  The choice of codons was calculated
% such that any combination of meaningful codons will never produce another
% meaningful codon in a different reading frame.
%
% This was disproved in 1961, when the first codon was discovered.  It
% mapped UUU (nonsense in a comma-free code) to the amino acid
% phenylalanine.  The majority of the genetic code had been mapped by 1965.
%
% The actual code ignores the issue of reading frames, instead focussing on
% error minimisation, which is achieved via a redundancy in the code.
% There is also a great deal of clustering, with similar codons usually
% coding for the same amino acid.  For example, alanine is coded for by any codon
% beginning 'GC', which allows the final base to mutate with no
% consequences.  The following function displays a chart of the specified
% genetic code, with each codon coloured according to the amino acid that
% it codes for.
%

figure;
gcView( 'Standard' );


%% Including the Cost of Mutations in the Similarity Measure.

[codons codonCount] = lookupTable( 'codons' );
[gcIndices gcCount] = lookupTable( 'gcIndices' );
names = {};
gcCosts = zeros( gcCount, 1 );
for gc = 1 : gcCount
    gcode = geneticcode( gcIndices(gc) );
    names{gc} = gcode.Name;
    for c = 1 : codonCount
        codon = codons{c};
        gcCosts(gc) = gcCosts(gc) + mutationCost( codon, gcode, 10 );
    end
end

fig = figure;
axes( 'Parent', fig, 'YTickLabel', names, 'YAxisLocation', 'right',...
     'YTick', 1:gcCount, 'YDir', 'reverse', 'Position', [.05 .10 .7 .85],...
     'XLim', [1600 1900] );
hold( 'all' );
barh( gcCosts );
title( 'The Cost of Mutation for all Genetic Codes' );
xlabel( 'Cost' );

% SPEED THIS UP
% [gcDist gcNames] = gcDistance( 'mutationCostPam10' );
% gcTree = seqlinkage( gcDist, 'UPGMA', gcNames );
% 
% roots = { 'Bacterial' };
% for rootID = 1 : length( roots )
%     rootName = roots( rootID );
%     root = getbyname( gcTree, rootName );
%     gcTree = reroot( gcTree, root );
%     plot( gcTree, 'type', 'square' );
%     title( rootName );
% end


%% Theories on the Origins of Life
%
% Our Last Universal Commmon Ancestor is thought to have existed somewhere
% in the region of three to four billion years ago [History of the Universe].
% It is thought that the Central Dogma of Molecular Biology, 'DNA makes RNA 
% makes Protein', applied even then.  This gives rise to the question of 
% what came first, DNA, RNA or proteins?
%
% Roughly speaking, DNA provides a reliable way to store information and
% proteins are repsonsible for metabolism, which constructs and maintains
% the cell, as well as producing energy.  RNA is not as good at long term
% information storage, as it lacks the complementary strand and
% also uses the base uracil, which is a lot more likely to mutate than
% thymine in DNA.
%
% This is a basis of the 'RNA World' hypothesis [Dworkin 2001] ...
%


%% Theories on the Evolution of the Genetic Code
%
% It is sensible to ask how and why the system of codons and translational
% apparatus came about.  A system with such complexity is extremely
% unlikely to have just appeared (unless resorting to intelligent design),
% so it must have developed from a simpler arrangement.  There are a number
% of theories as to how the genetic codes originated how they have since
% evolved.  
%
% The early theories were the 'Frozen Accident' [Crick 68] and
% stereochemical theories.  Stereochemical theories predict an affinity
% between particular amino acids and particular codons, and give a good
% explanation for the early appearance of such a code.  The frozen accident
% suggests that simpler codes were gradually altered in order to
% accommodate new amino acids and that the standard code has fixed due to
% the difficulty for it to change.
%
% More on the theories...
% 
% Trifonov

%% WEKA
%
% The following code outputs a WEKA format '.arff' file that contains data
% about each genetic code.  It is great for playing around with the data
% and looking for trends.  One problem was that it is set up to find
% general trends in the data, which would invariably ouput a copy of the
% standard genetic code.  So during use it is worth changing the setting
% of classifiers so that they overfit the data as far as possible.  This
% then finds the specific differences between genetic codes.  This is best
% visualised in the following section.
%

[codons codonCount] = lookupTable( 'codons' );
[gcIndices gcCount] = lookupTable( 'gcIndices' );

arffFile = 'geneticCodes.arff';
arff = fopen( arffFile, 'w' );
if arff == -1
    disp( [ 'Failed to open ' arffFile ]);
else
    fprintf( arff, '@relation geneticCodes\n\n' );
    fprintf( arff, '@attribute codon {' );
    for i = 1 : codonCount - 1
        fprintf( arff, '%s, ', char(codons(i)) );
    end
    fprintf( arff, '%s }\n', char(codons(codonCount)) );
    fprintf( arff, '@attribute codon1 {A, C, T, G}\n' );
    fprintf( arff, '@attribute codon2 {A, C, T, G}\n' );
    fprintf( arff, '@attribute codon3 {A, C, T, G}\n' );
    fprintf( arff, '@attribute protein { *, A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V }\n' );
    fprintf( arff, '@attribute codeName { "Standard", "Vertebrate Mitochondrial", "Yeast Mitochondrial", "Mold, Protozoan, and Coelenterate Mitochondrial and Mycoplasma/Spiroplasma", "Invertebrate Mitochondrial", "Ciliate, Dasycladacean and Hexamita Nuclear", "Echinoderm Mitochondrial", "Euplotid Nuclear", "Bacterial and Plant Plastid", "Alternative Yeast Nuclear", "Ascidian Mitochondrial", "Flatworm Mitochondrial", "Blepharisma Nuclear;", "Chlorophycean Mitochondrial", "Trematode Mitochondrial", "Scenedesmus obliquus mitochondrial", "Thraustochytrium Mitochondrial" }\n\n' );
    fprintf( arff, '@data\n\n' );
    
    for gc = 1 : gcCount
        code = geneticcode( gcIndices(gc) );
        for c = 1 : codonCount
            codon = codons{c};
            for i = 1 : 3
                bases(i) = codon(i);
            end
            protein = code.(codon);
            name = code.Name;
            fprintf( arff, '%s, %c, %c, %c, %c, "%s"\n', codon, bases(1), bases(2), bases(3), protein, name );
        end
    end
    fclose( arff );
    disp( [ 'Written genetic code data to ' arffFile ] );
end

%% Codon Changes
%
% This creates a figure showing the amount of variation in amino acids
% coded for by each codon, across all genetic codes.  It clearly shows a
% large amount of conservation in the codons beginning with G and C, (I
% think there is a reason, linked with other studies of GC content).
%

[codons codonCount] = lookupTable( 'codons' );
[gcIndices gcCount] = lookupTable( 'gcIndices' );
codonString = '';
for c = 1 : codonCount
    codon = codons{c};
    for gc = 1 : gcCount
        gcode = geneticcode( gcIndices(gc) );
        aminoAcids(gc) = gcode.(codon);        
    end
    % the variation is the number of unique proteins greater than 1
    variation = length( unique(aminoAcids) ) - 1;
    for i = 1 : variation
        codonString = [codonString codon];
    end
end
figure;
c = codoncount( codonString, 'frame', 1, 'figure', true );
title( 'Codon coding variation across genetic codes.' );


%% Ranked list of amino acids
% These are ranked by the number of codons that encode them.
% The reverse genetic code is looked up and the counts for each amino
% acid is added to the structure.
% This is then sorted and plotted.
%
aminos = {'Starts', 'Stops', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };
aminoCount = length( aminos );
gcNos = [ 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23 ]';
gcCount = length( gcNos );
ranked = zeros( aminoCount, gcCount );
for i = 1 : gcCount
    gc = revgeneticcode( gcNos(i) );
    % look up the length of the entry for each amino acid
    for aa = 1 : aminoCount
        aaField = aminos{aa};
        ranked(aa, i) = length( gc.(aaField) );
    end
end
% rank these entries
[y indices] = sort( ranked );
% plot
bar( y' )


%% Bibliography
%
% * Crick, Fancis. 1968. _The Origin of the Genetic Code_. Journal of
% Molecular Biology.
%
% % Dworkin, Jason P.  Lazcano, Antonio.  Miller, Stanley. 2001. _The roads
% to and from the RNA world_.  Journal of Theoretial Biology.
%
% * Gilis, Dimitri.  Massar, Serge.  Cerf, Nicolas J.  Rooman, Marianne.
% 2001. _Optimality of the genetic code with respect to protein stability 
% and amino-acid frequencies_.  Genome Biology.
%
% * Hayes, Brian. 1998. _The Invention of the Genetic Code_.  American
% Scientist.
%
% * History of the Universe - Timeline. http://www.historyoftheuniverse.com/tl1.html
%
%
%
%
%
%
%
%
%
%
%
%
%
