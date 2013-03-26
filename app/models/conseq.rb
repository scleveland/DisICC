class Conseq
  include DataMapper::Resource
  
  property :conseq_id, Serial
  property :seq_id, Integer, :required => true
  property :aasequence_id, Integer, :required => true
  property :position, Integer, :required => false
  property :score, Float, :required => true
  property :color, Integer, :required => false
  property :state, String, :required => false
  property :function, String, :required => false
  property :msa_data, String, :required => true
  property :residue_variety, String, :required => false
  property :deleted_at, ParanoidDateTime
  property :qq_interval, String, :required => false
  property :std, Float, :required => false
  
  #belongs_to :aasequence, 'AAsequence', :child_key => :aasequence_id
  #has 1, :sequence, :through=>:aasequence
end
