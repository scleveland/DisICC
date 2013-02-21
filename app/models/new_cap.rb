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
  property :greater_than_twenty_away, Boolean, :required => false
  property :greater_than_50_away, Boolean, :required => false
  property :greater_than_100_away, Boolean, :required => false
  property :greater_than_200_away, Boolean, :required => false
  property :greater_than_300_away, Boolean, :required => false
  property :greater_than_400_away, Boolean, :required => false
  property :greater_than_500_away, Boolean, :required => false
  property :greater_than_600_away, Boolean, :required => false
  property :greater_than_700_away, Boolean, :required => false
  property :greater_than_800_away, Boolean, :required => false
  property :greater_than_900_away, Boolean, :required => false
  property :greater_than_1000_away, Boolean, :required => false
  property :greater_than_1100_away, Boolean, :required => false
  property :greater_than_1200_away, Boolean, :required => false
  property :greater_than_1300_away, Boolean, :required => false
  property :greater_than_1400_away, Boolean, :required => false
  property :greater_than_1500_away, Boolean, :required => false
  property :greater_than_1600_away, Boolean, :required => false
  property :greater_than_1700_away, Boolean, :required => false
  property :greater_than_1800_away, Boolean, :required => false
  property :greater_than_1900_away, Boolean, :required => false
  property :greater_than_2000_away, Boolean, :required => false
  property :greater_than_2100_away, Boolean, :required => false
  property :greater_than_2200_away, Boolean, :required => false
  property :greater_than_2300_away, Boolean, :required => false
  property :greater_than_2400_away, Boolean, :required => false
  property :greater_than_2500_away, Boolean, :required => false
  # has n, :disorder_values
  # belongs_to :sequence
    belongs_to :sequence, 'Sequence', :child_key =>[:seq_id], :parent_key => [:seq_id]
    
  # eval_res_distances 
  # determine is the prediction's position one is further away than twenty residues from position two  
  def eval_res_distances
    self.greater_than_twenty_away = (self.position_one - self.position_two).abs > 20
    self.greater_than_50_away = (self.position_one - self.position_two).abs > 50
    self.greater_than_100_away = (self.position_one - self.position_two).abs > 100
    self.greater_than_200_away = (self.position_one - self.position_two).abs > 200
    self.greater_than_300_away = (self.position_one - self.position_two).abs > 300
    self.greater_than_400_away = (self.position_one - self.position_two).abs > 400
    self.greater_than_500_away = (self.position_one - self.position_two).abs > 500
    self.greater_than_600_away = (self.position_one - self.position_two).abs > 600
    self.greater_than_700_away = (self.position_one - self.position_two).abs > 700
    self.greater_than_800_away = (self.position_one - self.position_two).abs > 800
    self.greater_than_900_away = (self.position_one - self.position_two).abs > 900
    self.greater_than_1000_away = (self.position_one - self.position_two).abs > 1000
    self.greater_than_1100_away = (self.position_one - self.position_two).abs > 1100
    self.greater_than_1200_away = (self.position_one - self.position_two).abs > 1200
    self.greater_than_1300_away = (self.position_one - self.position_two).abs > 1300
    self.greater_than_1400_away = (self.position_one - self.position_two).abs > 1400
    self.greater_than_1500_away = (self.position_one - self.position_two).abs > 1500
    self.greater_than_1600_away = (self.position_one - self.position_two).abs > 1600
    self.greater_than_1700_away = (self.position_one - self.position_two).abs > 1700
    self.greater_than_1800_away = (self.position_one - self.position_two).abs > 1800
    self.greater_than_1900_away = (self.position_one - self.position_two).abs > 1900
    self.greater_than_2000_away = (self.position_one - self.position_two).abs > 2000
    self.greater_than_2100_away = (self.position_one - self.position_two).abs > 2100
    self.greater_than_2200_away = (self.position_one - self.position_two).abs > 2200
    self.greater_than_2300_away = (self.position_one - self.position_two).abs > 2300
    self.greater_than_2400_away = (self.position_one - self.position_two).abs > 2400
    self.greater_than_2500_away = (self.position_one - self.position_two).abs > 2500
    self.save
  end
  
end
