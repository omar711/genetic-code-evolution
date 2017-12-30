%% LOOKUPTABLE returns the specified table and its size.
% The following ids can be used:
%   'ntAlphabet'        - nucleotide alphabet
%   'aaAlphabet'        - amino acid alphabet
%   'codons'            - the 64 possible codons
%   'gcIndices'         - indices for each genetic code
function [info size] = lookupTable( id )
    switch lower( id )
        case 'ntalphabet'
            info = getNTalphabet;
        case 'aaalphabet'
            info = getAAalphabet;
        case 'codons'
            info = getCodons;
        case 'gcindices'
            info = getGcNos;
        otherwise
            error( [ 'Unrecognised lookup id: ' id ] );
    end
    size = length( info );

%%
% Get the genetic code numbers
function codes = getGcNos
    codes = [ 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23 ]';

%%
% Get the alphabet - array of characters.
function a = getNTalphabet
    a = ['A' 'G' 'C' 'T' ]';
    
%%
% Get the amino acid alphabet
function a = getAAalphabet
    a = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]';
    
%%
% Make a list of all possible codons
function codons = getCodons
    alphabet = getNTalphabet;
    aSize = length( alphabet );
    codons = { aSize * aSize * aSize };
    codon = 1;
    for i = 1 : aSize
        li = alphabet(i);
        for j = 1 : aSize
            lj = alphabet(j);
            for k = 1 : aSize
                lk = alphabet(k);
                codons{codon} = strcat( li, lj, lk );
                codon = codon + 1;
            end
        end
    end
    
