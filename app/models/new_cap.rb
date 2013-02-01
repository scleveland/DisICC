class NewCap
  include DataMapper::Resource
  
  property :caps_id, Serial
  property :aasequence_id, Integer, :required => true
  property :position_one, Integer, :required => true #this is the position of the aasequence
  property :position_two, Integer, :required => true # this is the position of the partner
  property :mean_one, Float, :required => true
  property :mean_two, Float, :required => true
  property :correlation, Float, :required => true
  property :seq_id, Integer, :required => true
  property :create_at, DateTime, :required => false
  # has n, :disorder_values
  # belongs_to :sequence
    belongs_to :sequence, 'Sequence', :child_key =>[:seq_id]
end
