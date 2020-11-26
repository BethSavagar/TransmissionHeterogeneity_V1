

effcIdGenerator <- function(effc, effc_id, N) {
  Vsize <- effc # vector Vsize contains the number of effective contacts for each infected unit(where effective contacts are non-zero)
  
  # Sample from all units in population to generated effectively contacted unit ids.
  # Use I_index_effc to remove the possibility of a self-self contact.
  # The list (L) contains, for each infected unit with non-zero effective contacts (effc/Vsize), the id of contacted units simulated by random sampling of the population.
  L <- lapply(1:length(Vsize),
              function(x)
                sample((1:N)[-effc_id[x]], # Sample contacted unit ids from the total population (N), excluding the current infected unit (defined by effc_id[x])
                       Vsize[x], # The number of ids sampled is defined by the number of effective contacts made by the 'x' Infected unit.
                       replace = FALSE)) # Once contacted a unit cannot be contacted again by the same infected unit.
  
  # Transform the list of contact ids into a vector which contains no duplicate ids (i.e. if a unit is contacted by two different infected units it is only recorded once, since a unit cannot be infected twice in 1 timestep)
  
  contact_ids <- unique(unlist(L)) # contact_ids is a vector which stored the ids of all units which are effectively contacted by an infected unit
  return(contact_ids)
}

  
  
  
