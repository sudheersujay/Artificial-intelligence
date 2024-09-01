myFunction = function(moveInfo, readings, positions, edges, probs) {
  currentNode = positions[3]

  bp1 = positions[1]
  if(is.na(bp1))
    bp1 = 0

  bp2 = positions[2]
  if(is.na(bp2))
    bp2 = 0

  st = rep(0, 40)
  e = getEmissionMatrix(probs, readings)

  if (moveInfo$mem$status < 2) {
    s0 = getInitialState(bp1, bp2, currentNode)
    moveInfo$mem$s0 = s0
    moveInfo$mem$status = 2

    for (i in 1:40) {
      #maybe add weight?
      st[i] = s0[i] * e[i]
    }
    s0 = st / sum(st) # normalize
    #print(st)

    crocLocation = which(st == max(st))
  }
  else if(bp1 < 0 || bp2 < 0){
    if (bp1 < 0) {
      crocLocation = -1 * bp1
      st[crocLocation] = 1
    }
    if (bp2 < 0) {
      crocLocation = -1 * bp2
      st[crocLocation] = 1
    }
  } else {
    s0 = moveInfo$mem$s0
    tMatrix = getTransitionMatrix(edges)
    st = hiddenMarkovModel(s0,e,tMatrix)
    crocLocation = which(st == max(st))
  }
  path = shortestBFS(currentNode, crocLocation, edges)
  if(length(path) >= 2) { # croc is two or more nodes away
    moveInfo$moves = c(path[1], path[2])
  } else if(length(path) == 1) { # croc is one node away
    moveInfo$moves = c(path[1], 0)
    st[crocLocation] = 0 # In case search misses, mark pool as not having croc
  } else if(length(path) == 0) { # croc is at same position
    moveInfo$moves=c(0,0)
    st[crocLocation] = 0 # In case search misses, mark pool as not having croc
  }

  moveInfo$mem$s0 = st

  return (moveInfo)
}

hiddenMarkovModel = function(stateMatrix,emissionMatrix,tMatrix)
{
  bestState = rep(0, 40)
  stateT = stateMatrix%*%tMatrix
  for(i in 1:40)
  {
    total = stateT[i] * emissionMatrix[i]
    bestState[i] = total
  }
  bestState = bestState/sum(bestState)
  return (bestState)
}

getInitialState = function(bp1, bp2, rangerPos, numNodes = 40) {
  # We assume backpackers and the ranger won't start on the same spot
  remNodes = numNodes-3
  s0 = rep(1/remNodes, numNodes)

  # Croc doesn't start on any of the backpackers or ranger's positions
  s0[rangerPos] = 0
  s0[bp1] = 0
  s0[bp2] = 0

  return (s0)
}

getTransitionMatrix = function(edges)
{
  tMatrix = matrix(, nrow = 40, ncol = 40)
  options(max.print=999999)
  for (i in 1:40)
  {
    for(j in 1:40)
    {
      op=c(edges[which(edges[,1]==i),2],edges[which(edges[,2]==i),1],i)
      len = length(op)
      if(j %in% op)
      {
        tMatrix[i,j] = 1/len
      }
      else
      {
        tMatrix[i,j] = 0
      }
    }
  }
  return(tMatrix)
}

getEmissionMatrix = function(probs, readings) {
  salinity = dnorm(readings[1],probs$salinity[,1],probs$salinity[,2])
  phosphate = dnorm(readings[2],probs$phosphate[,1],probs$phosphate[,2])
  nitrogen = dnorm(readings[3],probs$nitrogen[,1],probs$nitrogen[,2])
  emissionMatrix = salinity * phosphate * nitrogen
  emissionMatrix = emissionMatrix / sum(emissionMatrix) # normalize
  return (emissionMatrix)
}

shortestBFS = function(src, dest, edges)
{
  visited = c()
  frontier = c()
  parents = replicate(40, 0)

  currNode = src
  parents[src] = -1

  while(currNode != dest) {
    visited = append(visited, currNode)
    neighbours = getOptions(currNode, edges)
    for (node in neighbours) {
      if (!is.element(node, frontier) && !is.element(node, visited)) {
        frontier = append(frontier, node)
        parents[node] = currNode
      }
    }
    currNode = frontier[1]
    frontier = setdiff(frontier, currNode)
  }
  prevNode = dest
  nextNode = parents[dest]
  path = c()
  while(nextNode != -1) {
    path = c(c(prevNode), path)
    prevNode = nextNode
    nextNode = parents[nextNode]
  }

  return (path)
}
