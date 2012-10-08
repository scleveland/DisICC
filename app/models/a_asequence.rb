class AAsequence 
  include DataMapper::Resource
  
  property :AAsequence_id, Serial
  property :seq_id, Integer, :required => true
  property :amino_acid, String, :length=> 1,  :required => true
  property :original_position, Integer, :required => false
  property :disorder_consensus, Float, :required => false, :default => 0.0
  property :contact_consensus, Float, :required => false, :default => 0.0
  property :contact_positive_consensus, Integer, :required => false, :default => 0.0
  
  #belongs_to :sequence, 'Sequence', :parent_key=> [:seq_id], :child_key => [:seq_id]
  #has n, :disorder, 'Disorder', :parent_key=>[:disorder_id]
  #has n, :disorder_values, 'DisorderValue', :parent_key => disorder_value_id
end
