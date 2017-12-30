%% Return a measure of the cost of mutation for a given codon and geneticcode.
% This considers transitions and transversions seperately, assiging
% probabilities based on [Gilis 2001]:
% 
%   p(c' | c) = 1 / N       if they differ in 3rd base only.
%   p(c' | c) = 1 / N       for a transition in the 1st base.
%   p(c' | c) = 0.5 / N     for a transversion in the 1st base.
%   p(c' | c) = 0.5 / N     for a transition in the 2nd base.
%   p(c' | c) = 0.1 / N     for a transversion in the 2nd base.
%   p(c' | c) = 0           if more than 1 base differs.
%   where N is a normalisation factor to ensure that the 
%   sum of all p(c' | c) = 0, so N = 3.1
%
% The function enumerates each of these mutations and applies a cost
% measure from the PAM matrix.  Pass in the number of the pam matrix to use
% in the pamNum variable.
%
% Reference:
% * Gilis, Dimitri.  Massar, Serge.  Cerf, Nicolas J.  Rooman, Marianne.
% 2001. _Optimality of the genetic code with respect to protein stability 
% and amino-acid frequencies_.  Genome Biology.
%
function cost = mutationCost( codon, gcode, pamNum )
    N = 3.1;
    cost = 0;
    originalAA = gcode.(codon);
    
    % 3rd base differences and 1st base transitions
    prob = 1 / N;
    mutations = [ pointMutations( codon, 3 )...
        selectType( pointMutations( codon, 1 ), 'Transition' ) ];
    cost = cost + sumCost( mutations, prob, originalAA, gcode, pamNum );
    
    % 1st base transversions and 2nd base transitions
    prob = 0.5 / N;
    mutations = [ selectType( pointMutations( codon, 1 ), 'Transversion' )...
        selectType( pointMutations( codon, 2 ), 'Transition' ) ];
    cost = cost + sumCost( mutations, prob, originalAA, gcode, pamNum );
    
    % 2nd base transversions
    prob = 0.1 / N;
    mutations = selectType( pointMutations( codon, 2 ), 'Transversion' );
    cost = cost + sumCost( mutations, prob, originalAA, gcode, pamNum );
    
    
%% Sum the cost of mutations in the given list given:
% the probability of the mutation.
% the original amino acid.
% the genetic code.
% the number of the pam matrix for cost lookup.
%
function cost = sumCost( mutations, prob, originalAA, gcode, pamNum )
    cost = 0;
    for i = 1 : length( mutations )
        mutation = mutations{i};
        mutatedAA = gcode.( mutation.Codon );
        cost = cost + prob * ( pamLookup( pamNum, originalAA, mutatedAA ) );
    end
    
%% Return all single point mutations for a particular codon
% The field Type details the type of mutation and Codon gives the mutated
% codon.
function mutations = allMutations( codon )
    mutations = {};
    for cNt = 1 : length( codon )
        mutations = [ mutations pointMutations( codon, cNt ) ];
    end
    
%% Return the single point mutations at a particular point in a codon
% Also has a field to indicate whether it is a 'Transition' or a
% 'Transversion'
%
function pms = pointMutations( codon, point )
    [nts ntCount] = lookupTable( 'ntAlphabet' );

    index = 1;
    codonNt = codon(point);
    pms = {};
    for nt = 1 : ntCount
        mutantNt = nts(nt);
        if ( mutantNt ~= codonNt )
            mutated = codon;            
            mutated(point) = nts(nt);
            pms{index}.Codon = mutated;
            pms{index}.Type = mutationType( codonNt, mutantNt );
            index = index + 1;
        end
    end
    
%% Select the entries that are of the specified Type
function selected = selectType( mutations, type )
    % check for correct type input
    if ( strcmpi( type, 'transversion' ) ~= 0 && strcmpi( type, 'transition' ) ~= 0 )
        error( 'type should equal "transition" or "transversion"' );
    end
    
    selected = {};
    index = 1;
    for i = 1 : length( mutations )
        if ( strcmpi( mutations{i}.Type, type ) )
            selected{index} = mutations{i};
            index = index + 1;
        end
    end
        

%% Classify the type of mutation, 'Transition' or 'Transversion'
% nt2int returns odd for purine and even for pyrimidine, so this gets the
% nt2int for each nucleotide, if they match, 'Transition' is returned,
% otherwise 'Transversion' is returned.
function type = mutationType( ntOriginal, ntMutant )
    type = 'Transversion';
    original = mod( nt2int( ntOriginal ), 2 );
    mutant = mod( nt2int( ntMutant ), 2 );
    if ( original == mutant ), type = 'Transition'; end

    
    
%% Utility to lookup a distance in a pam matrix by letter.
% This also inverts the values in the pam matrix so that unrelated amino
% acids are large positive values and related ones are negative.
% Finally it adds an offset to make sure everything is positive.
function pamdist = pamLookup( pamNum, fromAA, toAA )
    [pamMatrix pamDetails] = pam( pamNum );
    order = pamDetails.Order;
    fromIndex = findstr( order, fromAA );
    toIndex = findstr( order, toAA );
    % invert and offset the value
    offset = max( max( pamMatrix ) );
    pamdist = pamMatrix( fromIndex, toIndex );
    pamdist = -pamdist + offset;
    
