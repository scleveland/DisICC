class IntraResidueContact
  include DataMapper::Resource
  
  property :intra_residue_contact_id, Serial
  property :seq_id, Integer, :required => true
  property :first_residue, Integer, :required => true
  property :second_residue, Integer, :required => true
  property :confidence, Float, :required => true
  property :type, String, :required => true

  

end
