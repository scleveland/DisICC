class Conseq
  include DataMapper::Resource
  
  property :conseq_id, Serial
  property :seq_id, Integer, :required => true
  property :aasequence_id, Integer, :required => true
  property :score, Float, :required => true
  property :color, Integer, :required => true
  property :state, String, :required => true
  property :function, String, :required => true
  property :msa_data, String, :required => true
  property :residue_variety, String, :required => false
  
  
end
