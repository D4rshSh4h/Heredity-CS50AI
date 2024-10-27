import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    set_people = set(people.keys())
    set_no_genes = set_people - one_gene - two_genes
    set_no_trait = set_people - have_trait
    prob_two_genes_with_trait = 1
    prob_two_genes_without_trait = 1
    prob_one_gene_with_trait = 1
    prob_one_gene_without_trait = 1
    prob_no_genes_with_trait = 1
    prob_no_genes_without_trait = 1

    
    for person in two_genes:
        mother = people[person]['mother']
        father = people[person]['father']
        if mother == None and father == None:
            probability = PROBS["gene"][2]
        else:
            if mother in two_genes and father in two_genes:
                probability = 2*(1-PROBS["mutation"])
            elif mother in set_no_genes and father in set_no_genes:
                probability = PROBS["mutation"]**2
            elif mother in one_gene and father in one_gene:
                probability = 0.5*2*(1-PROBS["mutation"])
            elif mother in one_gene and father in set_no_genes:
                probability = 0.5*(1-PROBS["mutation"])*PROBS["mutation"]
            elif father in one_gene and mother in set_no_genes:
                probability = 0.5*(1-PROBS["mutation"])*PROBS["mutation"]
            elif mother in two_genes and father in set_no_genes:
                probability = (1-PROBS["mutation"])*PROBS["mutation"]
            elif father in two_genes and mother in set_no_genes:
                probability = (1-PROBS["mutation"])*PROBS["mutation"]
            elif mother in one_gene and father in two_genes:
                probability = 0.5*(1-PROBS["mutation"])**2
            elif father in one_gene and mother in two_genes:
                probability = 0.5*(1-PROBS["mutation"])**2
            else:
                raise TypeError

                
        if person in have_trait:
            probability *= PROBS["trait"][2][True]
            prob_two_genes_with_trait *= probability
        elif person in set_no_trait:
            probability *= PROBS["trait"][2][False]
            prob_two_genes_without_trait *= probability
        else:
            print('Error')

    
    for person in one_gene:
        mother = people[person]['mother']
        father = people[person]['father']
        if mother == None and father == None:
            probability = PROBS["gene"][1]
        else:
            if mother in two_genes and father in two_genes:
                probability = 1
            elif mother in set_no_genes and father in set_no_genes:
                probability = PROBS["mutation"]*(1-PROBS["mutation"])*2
            elif mother in one_gene and father in set_no_genes:
                probability = (0.5*(1-PROBS["mutation"])*(1-PROBS["mutation"])) + (0.5*(1-PROBS["mutation"])*PROBS["mutation"])
            elif father in one_gene and mother in set_no_genes:
                probability = (0.5*(1-PROBS["mutation"])*(1-PROBS["mutation"])) + (0.5*(1-PROBS["mutation"])*PROBS["mutation"])
            elif mother in two_genes and father in set_no_genes:
                probability = ((1-PROBS["mutation"])**2)+(PROBS["mutation"])**2
            elif father in two_genes and mother in set_no_genes:
                probability = ((1-PROBS["mutation"])**2)+(PROBS["mutation"])**2
            elif mother in one_gene and father in one_gene:
                probability = ((0.5*(1-PROBS["mutation"]))**2)*2
            elif mother in one_gene and father in two_genes:
                probability = ((1-PROBS["mutation"])**2)*0.5*2
            elif mother in two_genes and father in one_gene:
                probability = ((1-PROBS["mutation"])**2)*0.5*2
            else:
                raise TypeError
            
        if person in have_trait:
            probability *= PROBS["trait"][1][True]
            prob_one_gene_with_trait *= probability
        elif person in set_no_trait:
            probability *= PROBS["trait"][1][False]
            prob_one_gene_without_trait *= probability
        else:
            print('Error')
    
    for person in set_no_genes:
        mother = people[person]['mother']
        father = people[person]['father']
        if mother == None and father == None:
            probability = PROBS["gene"][0]
        else:
            if mother in two_genes and father in two_genes:
                probability = (PROBS["mutation"])**2
            elif mother in set_no_genes and father in set_no_genes:
                probability = (1-PROBS["mutation"])**2
            elif mother in one_gene and father in set_no_genes:
                probability = 0.5*(1-PROBS["mutation"])
            elif father in one_gene and mother in set_no_genes:
                probability = 0.5*(1-PROBS["mutation"])
            elif mother in two_genes and father in set_no_genes:
                probability = (PROBS["mutation"])*(1-PROBS["mutation"])
            elif father in two_genes and mother in set_no_genes:
                probability = (PROBS["mutation"])*(1-PROBS["mutation"])
            elif mother in one_gene and father in one_gene:
                probability = (0.5*(1-PROBS["mutation"]))**2
            elif mother in one_gene and father in two_genes:
                probability = (PROBS["mutation"])*0.5*(1-PROBS["mutation"])
            elif mother in two_genes and father in one_gene:
                probability = (PROBS["mutation"])*0.5*(1-PROBS["mutation"])
            else:
                raise TypeError

        
        if person in have_trait:
            probability *= PROBS["trait"][0][True]
            prob_no_genes_with_trait *= probability
        elif person in set_no_trait:
            probability *= PROBS["trait"][0][False]
            prob_no_genes_without_trait *= probability
        else:
            print('Error')
    
    return prob_no_genes_without_trait*prob_no_genes_with_trait*prob_one_gene_with_trait*prob_one_gene_without_trait*prob_two_genes_with_trait*prob_two_genes_without_trait


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        if person in one_gene:
            if person in have_trait:
                probabilities[person]["gene"][1] += p
                probabilities[person]["trait"][True] += p
            else:
                probabilities[person]["gene"][1] += p
                probabilities[person]["trait"][False] += p
        elif person in two_genes:
            if person in have_trait:
                probabilities[person]["gene"][2] += p
                probabilities[person]["trait"][True] += p
            else:
                probabilities[person]["trait"][False] += p
                probabilities[person]["gene"][2] += p
        elif person not in one_gene and person not in two_genes:
            if person in have_trait:
                probabilities[person]["gene"][0] += p
                probabilities[person]["trait"][True] += p
            else:
                probabilities[person]["gene"][0] += p
                probabilities[person]["trait"][False] += p
        else:
            if person in have_trait:
                probabilities[person]["trait"][True] += p
            else:
                probabilities[person]["trait"][False] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        for field in probabilities[person]:
            total = sum(probabilities[person][field].values())
            for value in probabilities[person][field]:
                probabilities[person][field][value] = ((probabilities[person][field][value])*1)/total


if __name__ == "__main__":
    main()
