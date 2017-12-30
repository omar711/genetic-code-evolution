%% Coding Differences
% This prints out the differences between the codes, specified by
% aTitle and bTitle.
% Whenever a codon codes for a different amino acid, it adds a field to
% the returned structure with the name of the codon.  This contains the
% amino acid coded by code A followed by the one coded by code B.
% e.g.
%   ATA: {'I'  'M'} is a difference between the standard and the yeast
%   mitochondrial genetic codes.
% It also adds entries for both start codons, StartsA and StartsB.
function differences = codingDifferences( aTitle, bTitle )
    codeA = geneticcode( aTitle );
    codeB = geneticcode( bTitle );
    [codons codonCount] = lookupTable( 'codons' );
    differences = struct();
    for c = 1 : codonCount
        codon = codons{c};
        aaA = codeA.(codon); aaB = codeB.(codon);
        if aaA ~= aaB
            differences.(codon) = { aaA, aaB };
        end
    end
    
    differences.StartsA = codeA.Starts;
    differences.StartsB = codeB.Starts;
    
    disp( ['Coding differences between ' aTitle ' and ' bTitle ' genetic codes.'] );
    disp( differences );