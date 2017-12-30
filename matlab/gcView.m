%% Produce a colour map of a genetic code.
% This is similar to the method codoncount, except that codons are coloured
% according to the amino acid that they code for.
%
% It may be worth adding options for the colour order, e.g. to show the
% hydrophobicity of amino acids.
%
function gcView( codeName )
    % assign colours to each codon
    [codons codonCount] = lookupTable( 'codons' );
    gCode = geneticcode( codeName );
    codonMap = zeros( 8, 8 );
    for c = 1 : codonCount
        codon = codons{c};
        [i j] = calculateIndex( c );
        codonMap(i, j) = aa2int( gCode.(codon) );
    end
    
    % dim the hsv colour map
    colours = jet * 0.9;
    
    % draw the figure
    imagesc( codonMap );
    axis off;
    colormap( colours );
    title( ['Codon Map of the ' gCode.Name ' Genetic Code'] );
    for c = 1 : codonCount
        codon = codons{c};
        [i j] = calculateIndex( c );
        text( j, i, codon, 'color', 'w', 'horizontalAlignment', 'center' );
    end
    
%% Calculate an 2d index for the grid, given 1d index.
% puts them in blocks similar to those in codoncount except 
% that a's are placed near to g's.
function [i j] = calculateIndex( ind )
    % four boxes, 1, 2; 3, 4
    box = ceil( ind / 16 );
    boxJ = mod( box - 1, 2 ) + 1;
    boxI = ceil( box / 2 );
    % shift the index to 1:16
    internalIndex = ind - ( (box - 1) * 16 );
    % work out indices as though this was in a 4x4 box
    internalI = ceil( internalIndex / 4 );
    internalJ = mod( internalIndex - 1, 4 ) + 1;
    % now offset these based on the box
    i = internalI + ( (boxI - 1) * 4 );
    j = internalJ + ( (boxJ - 1) * 4 );
    