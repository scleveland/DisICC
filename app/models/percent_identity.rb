class PercentIdentity
  include DataMapper::Resource
  
  property :seq1_id,  Integer, :key => true
  property :seq2_id, Integer, :key => true
  property :alignment_name, String, :key => true
  property :percent_id, Float, :required=> true
  property :deleted_at, ParanoidDateTime
  
  def seq2_sequence
    Sequence.get(seq2_id)
  end
end