class AlignmentPosition
  include DataMapper::Resource
  
  property :alignment_position_id, Serial
  property :alignment_id, Integer, :required => true
  property :position, Integer, :required => true
  property :aasequence_id, Integer, :required => true
  property :deleted_at, ParanoidDateTime
  
  # has n, :disorder_values
  # belongs_to :sequence
  
end
