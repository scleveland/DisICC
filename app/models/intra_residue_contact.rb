class IntraResidueContact
  include DataMapper::Resource
  
  property :id, Serial, :field=>'intra_residue_contact_id'
  property :seq_id, Integer, :required => true
  property :first_residue, Integer, :required => true
  property :second_residue, Integer, :required => true
  property :confidence, Float, :required => true
  property :type, String, :required => true

  
  belongs_to :sequence, 'Sequence', :child_key => :seq_id
end
