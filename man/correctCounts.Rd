\name{correctCounts}
\alias{correctCounts}
\title{
  Correcting multiply mapped counts
}
\description{
using a base count file (consist of columns: id, count) and a multiply-mapped fragment file (consist of columns: fragment_id, gene_id )

}
\usage{
correctCounts(baseCount, multiCount)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseCount}{
    basecount should be a dataframe for storing gene count from uniquely-mapped reads
    Have to contain columns named as 'id' for gene_id
      and 'count' for counts
  }
  \item{multiCount}{
    multiCount should be a dataframe for storing multiply-mapped reads and their mapped locus.
    columns named as 'fragment_id' for storing read id
      and 'gene_id' for storing gene_id
  }
}

\author{
Douglas wu
}

\examples{

library(correctMultiCount)

data(baseCount)
data(multiCount)

## The function is currently defined as
df = correctCounts(baseCount,multiCount)
head(df)
}
